# %%
import io
import json
import subprocess
from functools import partial, reduce
from itertools import cycle, permutations
from pathlib import Path
from typing import Iterable

import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
from expression import pipe

from oligocheck.merfish.alignment import run_bowtie
from oligocheck.merfish.external_data import ExternalData
from oligocheck.merfish.filtration import handle_overlap
from oligocheck.merfish.readouts.blacklist import check_align, get_blacklist
from oligocheck.merfish.readouts.codebook import CodebookPicker

# from oligocheck.merfish.filtration import the_filter
from oligocheck.sequtils import equal_distance, parse_cigar, parse_sam, reverse_complement, tm_hybrid

pl.Config.set_fmt_str_lengths(30)
pl.Config.set_tbl_rows(30)
wants = ["Eef2"]


def count_genes(df: pl.DataFrame) -> pl.DataFrame:
    return df.groupby("gene").count().sort("count")


def filter_gene(df: pl.DataFrame, gene: str) -> pl.DataFrame:
    return df.filter(pl.col("gene") == gene)


# def handle_done(df: p.DataFrame, target: int = 40):
#     if len(df) < target:
#         return df
#     picked = df.sort_values(by="pos", ascending=False)
#     return picked.iloc[equal_distance(len(picked), target)]


ensembl = ExternalData(
    "data/mm39/gencode_vM32_transcripts_all.parquet",
    path="data/mm39/Mus_musculus.GRCm39.109.gtf",
    fasta="data/mm39/combi.fa",
)

# %%
temp = []
for gene in wants:
    try:
        temp.append(pl.read_parquet(f"output/{gene}.parquet"))
    except FileNotFoundError:
        print("File not found", gene)

metadata = [json.loads(Path(f"output/{gene}.json").read_text()) for gene in wants]
allowed = {x["gene"]: x["allowed"] for x in metadata}
df = pl.concat(temp).with_columns(pl.sum(pl.col("^ok_.*$")).alias("oks"))


def the_filter(df: pl.DataFrame, overlap: int = -1):
    return df.groupby("gene").apply(
        lambda group: handle_overlap(
            ensembl,
            group,
            criteria=[
                # fmt: off
                    pl.col("tm").is_between(49, 54) & (pl.col("oks") > 4) & (pl.col("hp") < 35) & (pl.col("nonspecific_binding") < 0.001),
                    pl.col("tm").is_between(49, 54) & (pl.col("oks") > 4) & (pl.col("hp") < 40) & (pl.col("nonspecific_binding") < 0.05),
                    pl.col("tm").is_between(47, 56) & (pl.col("oks") > 3) & (pl.col("hp") < 40),
                    pl.col("tm").is_between(46, 56) & (pl.col("oks") > 3) & (pl.col("hp") < 40),
                    pl.col("tm").is_between(46, 56) & (pl.col("oks") > 2) & (pl.col("hp") < 40),
                # fmt: on
            ],
            overlap=overlap,
        )
    )


# %%
out_nool = the_filter(df)
# %%
counts = count_genes(out_nool)
easy = counts.filter(pl.col("count") >= 45)["gene"]
noteasy = counts.filter(pl.col("count") < 45)["gene"]
res = the_filter(df.filter(pl.col("gene").is_in(noteasy)), overlap=15)
out = pl.concat([out_nool.filter(pl.col("gene").is_in(easy)), res])
# %%
# wants_filtered = [x for x in wants.iloc[:110] if (x != "Stmn1") and x not in noteasy[-9:]]
# out = out[out.gene.isin(wants_filtered)]

# %%
fpkm = pl.read_parquet("data/fpkm/P0_combi.parquet")
fp = {
    name: value
    for name, value in (
        ensembl.filter(pl.col("gene_name").is_in(wants))
        .join(fpkm, left_on="transcript_id", right_on="transcript_id(s)")
        .groupby("gene_name")
        .agg(pl.sum("FPKM"))
        .sort(by="FPKM", descending=False)
        .iter_rows()
    )
}
fpfil = [fp[x] for x in wants]


# %%
mhd4 = CodebookPicker(np.loadtxt("./mhd4_16bit.csv", delimiter=",", dtype=bool))
best, values = mhd4.find_optimalish(fpfil)
print(best, values)

readouts = pd.read_csv("zeroth_readouts.csv", index_col=0)
blacklist = get_blacklist("zeroth_readouts.csv", ensembl)


# %%


def gen_full_probe(
    df: pd.DataFrame,
    readouts: pd.Series,
    genes: Iterable[str],
    mhd4: npt.NDArray[np.bool_],
    header: str,
    footer: str,
    seed: int = 0,
):
    df = df[df.gene.isin(genes)].copy()
    assert mhd4.shape[0] > mhd4.shape[1]
    for gene, codeset in zip(genes, mhd4):
        rand = np.random.RandomState(seed)
        probes = df[df.gene == gene]
        if len(probes) < 45:
            print("Low probes for ", gene, len(probes))
        code_ids = np.argwhere(codeset).squeeze()
        codeseqs = readouts.iloc[code_ids].to_numpy().squeeze()
        combi = []
        for cs in permutations(zip(code_ids, codeseqs), 4):
            cids = [x[0] for x in cs]
            if (cids[0], cids[1]) in blacklist or (cids[2], cids[3]) in blacklist:
                continue
            combi.append(cs)
        assert len(combi) > 4

        rand.shuffle(combi)
        for (name, probe), codes in zip(probes.iterrows(), cycle(combi)):
            df.loc[name, "constructed"] = reverse_complement(
                codes[0][1] + "TT" + codes[1][1] + "TT" + probe.seq + "TT" + codes[2][1] + "TT" + codes[3][1]
            )
            df.loc[name, ["code1", "code2", "code3", "code4"]] = [x[0] for x in codes]

    df.dropna(inplace=True)
    df["constructed"] = df["constructed"].apply(lambda x: header + x + footer)
    return df


# %%

header = "CAACCGTACCGCTTGCTTACC"  # ps_ir73
footer = "AGGCGGTT"  # ps_if15
# header = "GAGAGGCGAGGACACCTACAG"  # ps_or2

# %%
out_nool = the_filter(df).filter(pl.col("priority") < 4).sort("pos_start")


# %%
zeroth_readout = pl.read_csv("zeroth_readouts.csv")
ros = zeroth_readout[22][0, "seq"], zeroth_readout[23][0, "seq"]
base = lambda r, seq: seq + "AA" + r + "AA" + r  # noqa: E731
ross = [partial(base, r) for r in ros]


for_bowtie = (
    out_nool.with_row_count("id")
    .with_columns(pl.col("id").apply(lambda i: i // (len(out_nool) / 2)))
    .with_columns(
        [
            pl.when(pl.col("id") == 0)
            .then(header + pl.col("seq").apply(ross[0]) + "TATTTCCC")
            .otherwise(header + pl.col("seq").apply(ross[1]) + "TATTTCCC")
            .alias("constructed"),
        ]
    )
)

f = io.StringIO()
for row in for_bowtie.iter_rows(named=True):
    f.write(f"@{row['name']}\n")
    f.write(row["constructed"] + "\n")
    f.write("+\n")
    f.write("~" * len(row["constructed"]) + "\n")

sam = run_bowtie(f.getvalue(), "data/mm39/mm39")

# %%

# %%

# from 18 for cellcycle
for_bowtie = gen_full_probe(out, zeroth_readout, wants_filtered, y, header, "TATTTCCC")
# %%


# %%
sam = parse_sam(sam)

# %%


def parse_cigar(df: pl.DataFrame):
    """
    Extract the position and length of matches.
    Most common is SMS, does not deal with indels accurately.
    But if there are indels, that means that the match is good enough to be discarded if non-specific.
    """
    return df.with_columns(
        match_length=(
            pl.col("cigar")
            .str.extract_all(r"(\d+)M")
            .arr.eval(pl.element().str.extract(r"(\d+)").cast(pl.UInt16))
            .arr.max()
        ),
        match_start=pl.col("cigar").str.extract(r"(\d+)S").cast(pl.UInt16),
    ).with_columns(
        seq_aligned=(
            pl.struct(["seq", "match_length", "match_start"]).apply(
                lambda t: t["seq"][t["match_start"] : t["match_start"] + t["match_length"]]
            )
        )  # type: ignore
    )


def check_align(df: pl.DataFrame, allowed_ts: dict[str, list[str]], match_dump_threshold: int = 24):
    # check = df.filter(pl.col("match_length") > match_dump_threshold)

    # check.groupby("gene").agg(pl.apply(exprs=["gene", "seq_aligned"], f=))

    def calc_offtarget(group: pl.DataFrame):
        gene = group[0, "gene"]
        return (
            group.with_columns(
                tm_hybrid=pl.when(pl.col("transcript").is_in(allowed_ts[gene]))
                .then(pl.lit(0))
                .otherwise(pl.col("seq_aligned").apply(tm_hybrid))
            )
            .filter(pl.col("tm_hybrid") > 0)
            .groupby("name")
            .agg(pl.all().sort_by("tm_hybrid").last())
        )
        # ).filter(pl.col("tm_hybrid") == pl.col("tm_hybrid").max().over("name"))

    # assert worst_case.iloc[40].tm < 40
    return df.groupby("gene").apply(calc_offtarget)


checked = (
    check_align(parse_cigar(sam), allowed, 17)
    .with_columns(mapped_to=pl.col("transcript").apply(ensembl.ts_to_gene))
    .with_columns(
        bad_tm=pl.col("tm_hybrid") > 33,
        horrible_tm=pl.col("tm_hybrid") > 37,
    )
)

checked = (
    for_bowtie.join(
        checked.select(["name", "tm_hybrid", "mapped_to", "bad_tm", "horrible_tm"]), on="name", how="left"
    )
    .with_columns(
        pl.col("tm_hybrid").fill_null(0),
        pl.col("bad_tm").fill_null(False),
        pl.col("horrible_tm").fill_null(False),
    )
    .filter(~pl.col("horrible_tm"))
)
# %%


def final_filter(df: pl.DataFrame):
    sortd = df.sort(
        [
            "bad_tm",
            "priority",
            "n_mapped",
        ],
        descending=False,
    ).with_columns(
        final=pl.col("constructed") + "TATAGTGAGTCGTATTA" + footer
    )  # ps_ir29
    # sortd["constructed"] += "TATAGTGAGTCGTATTAGACCGGTCT"  # ps_ir49
    return sortd


finale = checked.groupby("gene").apply(final_filter)
# finale = finale.astype(dict(code1=int, code2=int, code3=int, code4=int))
# finale.to_csv("tricycle_probes.csv")

# %%
t7 = reverse_complement("TAATACGACTCACTATAGGG")

assert all(finale.constructed.str.find(t7) != -1)

codes = set()
for s1, row in finale.groupby("gene")[["code1", "code2", "code3", "code4"]].agg(["unique"]).iterrows():
    assert len(row) == 4
    gene_code = tuple(sorted(reduce(lambda x, y: x & y, map(set, row.values))))
    codes.add(gene_code)

assert len(codes) == finale.groupby("gene").count().shape[0]


# %%
count_genes(finale)
# %%
# %%# %%
