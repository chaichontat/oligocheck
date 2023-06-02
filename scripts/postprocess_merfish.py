# %%
import json
from functools import partial, reduce
from itertools import cycle, permutations
from pathlib import Path
from typing import Any, Callable, Collection, Iterable, Sequence

import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
from expression import pipe

from oligocheck.algorithms import find_overlap, find_overlap_weighted
from oligocheck.merfish.alignment import gen_fastq, run_bowtie
from oligocheck.merfish.external_data import ExternalData
from oligocheck.merfish.filtration import handle_overlap
from oligocheck.merfish.readouts.blacklist import get_blacklist
from oligocheck.merfish.readouts.codebook import CodebookPicker
from oligocheck.seqcalc import tm_q5

# from oligocheck.merfish.filtration import the_filter
from oligocheck.sequtils import equal_distance, parse_sam, reverse_complement, tm_hybrid

pl.Config.set_fmt_str_lengths(30)
pl.Config.set_tbl_rows(30)


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
    fasta="data/mm39/combi.fa.gz",
)

# %%
temp = []
wants = list(filter(lambda x: x, Path("celltypegenes.csv").read_text().splitlines()))
# wants = ["Eef2"]
for gene in wants:
    try:
        temp.append(pl.read_parquet(f"output/{gene}.parquet"))
    except FileNotFoundError:
        print("File not found", gene)

metadata = [json.loads(Path(f"output/{gene}.json").read_text()) for gene in wants]
allowed = {x["gene"]: x["allowed"] for x in metadata}
df = pl.concat(temp).with_columns(pl.sum(pl.col("^ok_.*$")).alias("oks"))

# %%


def handle_overlap(
    ensembl: ExternalData,
    df: pl.DataFrame,
    criteria: list[pl.Expr],
    overlap: int = -1,
    n: int = 200,
):
    if len(gene := df.select(pl.col("gene").unique())) > 1:
        raise ValueError("More than one gene in filtered")
    gene = gene.item()
    df = df.sort(by=["is_ori_seq", "transcript_ori", "pos_end", "tm"], descending=[True, False, False, True])
    eid = ensembl.gene_to_eid(gene)
    tss = tuple(ensembl.filter(pl.col("gene_id") == eid)["transcript_id"])

    if not criteria:
        criteria = [pl.col("*")]

    ddf = df.lazy().with_row_count("index").with_columns(priority=pl.lit(0, dtype=pl.UInt8))
    priority = len(criteria)
    for criterion in reversed(criteria):
        ddf = ddf.update(
            ddf.filter(criterion).with_columns(priority=pl.lit(priority, dtype=pl.UInt8)), on="index"
        )
        priority -= 1
    df = ddf.filter(pl.col("priority") > 0).collect()

    selected_global = set()
    for ts in tss:
        this_transcript = df.filter(pl.col("transcript_ori") == ts)
        if not len(this_transcript):
            continue

        for i in range(1, len(criteria) + 1):
            run = (
                df.filter((pl.col("priority") <= i) & ~pl.col("index").is_in(selected_global))
                .select(["index", "pos_start", "pos_end", "priority", "n_mapped"])
                .sort(["pos_end", "pos_start"])
            )

            priorities = (run["priority"].apply(lambda x: 5 - x) + run["n_mapped"] / 5).to_list()
            if i == 1:
                ols = find_overlap(run["pos_start"], run["pos_end"], overlap=overlap)
            else:
                ols = find_overlap_weighted(run["pos_start"], run["pos_end"], priorities, overlap=overlap)
            sel_local = set(run[ols]["index"].to_list())
            print(i, len(sel_local))

            # for idx, st, end in (
            #     df.filter((pl.col("priority") <= i) & ~pl.col("index").is_in(selected_global))
            #     .select(["index", "pos_start", "pos_end"])
            #     .sort(["pos_end", "pos_start"])
            #     .iter_rows()
            # ):
            #     if st > curr_right - overlap:
            #         curr_right = end
            #         sel_local.add(idx)
            if len(sel_local) > n:
                break

        selected_global |= sel_local

    return df.filter(pl.col("index").is_in(selected_global))


def the_filter(df: pl.DataFrame, overlap: int = -1):
    out = []
    for name, group in df.groupby("gene"):
        out.append(
            handle_overlap(
                ensembl,
                group,
                criteria=[
                    # fmt: off
                    pl.col("tm").is_between(49, 54) & (pl.col("oks") > 4) & (pl.col("hp") < 35) & (pl.col("nonspecific_binding") < 0.001),
                    pl.col("tm").is_between(49, 54) & (pl.col("oks") > 4) & (pl.col("hp") < 40) & (pl.col("nonspecific_binding") < 0.05),
                    pl.col("tm").is_between(48, 56) & (pl.col("oks") > 3) & (pl.col("hp") < 40),
                    pl.col("tm").is_between(48, 56) & (pl.col("oks") > 2) & (pl.col("hp") < 40),
                    # fmt: on
                ],
                overlap=overlap,
            )
        )
    return pl.concat(out)


# %%
out_nool = the_filter(df.filter(pl.col("gene") == "Eef2"))
# %%


def stripplot(*args: Iterable[float], **kwargs: Any):
    import pandas as pd
    import seaborn as sns

    sns.set()

    df = pd.concat(pd.DataFrame({"x": u, "y": str(i)}) for i, u in enumerate(args))
    return sns.stripplot(data=df, x="x", y="y", **(dict(orient="h", alpha=0.6) | kwargs))


counts = count_genes(out_nool)
easy = counts.filter(pl.col("count") >= 45)["gene"]
noteasy = counts.filter(pl.col("count") < 45)["gene"]
print(counts)

# %%
stripplot(out_nool["pos_start"], df.filter(pl.col("gene") == "Eef2")["pos_start"])
# %%
res = the_filter(df.filter(pl.col("gene").is_in(noteasy)), overlap=15)
print(count_genes(res))
# %%
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


def stitch(seq: str, codes: Sequence[str]) -> str:
    return codes[0] + "TT" + codes[1] + "TT" + seq + "TT" + codes[2] + "TT" + codes[3]


# def gen_full_probe(
#     df: pd.DataFrame,
#     readouts: pd.Series,
#     genes: Iterable[str],
#     mhd4: npt.NDArray[np.bool_],
#     header: str,
#     footer: str,
#     seed: int = 0,
#     f_stitch: Callable[[str, Sequence[str]], str] = stitch,
# ):
#     df = df[df.gene.isin(genes)].copy()
#     assert mhd4.shape[0] > mhd4.shape[1]
#     for gene, codeset in zip(genes, mhd4):
#         rand = np.random.RandomState(seed)
#         probes = df[df.gene == gene]
#         if len(probes) < 45:
#             print("Low probes for ", gene, len(probes))
#         code_ids = np.argwhere(codeset).squeeze()
#         codeseqs = readouts.iloc[code_ids].to_numpy().squeeze()
#         combi = []
#         for cs in permutations(zip(code_ids, codeseqs), 4):
#             cids = [x[0] for x in cs]
#             if (cids[0], cids[1]) in blacklist or (cids[2], cids[3]) in blacklist:
#                 continue
#             combi.append(cs)
#         assert len(combi) > 4

#         rand.shuffle(combi)
#         for (name, probe), codes in zip(probes.iterrows(), cycle(combi)):
#             # Reverse complemented here.
#             # Always use mRNA sequence and readouts with Cs.
#             df.loc[name, "constructed"] = reverse_complement(f_stitch(probe["seq"], [x[1] for x in codes]))
#             df.loc[name, ["code1", "code2", "code3", "code4"]] = [x[0] for x in codes]

#     df.dropna(inplace=True)
#     df["constructed"] = df["constructed"].apply(lambda x: header + x + footer)
#     return df


# %%

header = "CAACCGTACCGCTTGCTTACC"  # ps_ir73
footer = "AGGCGGTT"  # ps_if15
# header = "GAGAGGCGAGGACACCTACAG"  # ps_or2

# %%
out_nool = the_filter(df).filter(pl.col("priority") < 4).sort("pos_start")


# %%
zeroth_readout = pl.read_csv("zeroth_readouts.csv")
ros = zeroth_readout[22][0, "seq"], zeroth_readout[23][0, "seq"]
base = lambda r, seq: reverse_complement(r + "AT" + seq + "TA" + r)  # noqa: E731
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

fastq = gen_fastq(for_bowtie.select(name=pl.col("name"), seq=pl.col("constructed")))
sam = run_bowtie(fastq.getvalue(), "data/mm39/mm39")

# %%

# %%

# from 18 for cellcycle
# for_bowtie = gen_full_probe(out, zeroth_readout, wants_filtered, y, header, "TATTTCCC")
# %%


# %%


def check_specificity(
    s: pl.DataFrame | Collection[str],
    seed_length: int = 13,
    reference: str = "data/mm39/mm39",
    threshold: int = 15,
    n_return: int = -1,
):
    if isinstance(s, Collection):
        fastq = pl.DataFrame(dict(name=[f"{i:04d}" for i in range(len(s))], seq=s))
    else:
        fastq = s

    fastq = fastq.with_columns(pl.col("seq").str.replace_all(" ", ""))

    sam = run_bowtie(
        gen_fastq(fastq).getvalue(),
        reference,
        seed_length=seed_length,
        n_return=n_return,
        threshold=threshold,
    )
    t = parse_sam(sam, split_name=False)

    # if parse_all_cigar:
    #     t = t.join(
    #         t[["id", "cigar"]]
    #         .with_columns(cigar=pl.col("cigar").str.extract_all(r"(\d+)M"))
    #         .explode("cigar")
    #         .with_columns(pl.col("cigar").str.extract(r"(\d+)").cast(pl.UInt8))
    #         .groupby("id")
    #         .agg(match=pl.col("cigar").sum()),
    #         on="id",
    #         how="left",
    #     )
    # else:
    #     t = t.with_columns(match=pl.col("cigar").str.extract(r"(\d+)M").cast(pl.UInt8))

    y = (
        t.with_columns(
            transcript_id=pl.col("transcript").str.extract(r"(ENSMUST.+)\.")
            # .arr.eval(pl.element().cast(pl.UInt8))
            # .max()
        ).join(fpkm, left_on="transcript_id", right_on="transcript_id(s)", how="left")
        # .with_columns(matchname=pl.col("transcript").apply(gtf_all.ts_to_gene), tm=pl.col("seq").apply(tm_fish))
    )

    y = y.join(
        y[["id", "mismatched_reference"]]
        .with_columns(mismatched_reference=pl.col("mismatched_reference").str.extract_all(r"(\d+)"))
        .explode("mismatched_reference")
        .with_columns(pl.col("mismatched_reference").cast(pl.UInt8))
        .groupby("id")
        .agg(
            match=pl.col("mismatched_reference").sum(),
            match_max=pl.col("mismatched_reference").max(),
        ),
        on="id",
        how="left",
    )

    y = y.with_columns(
        worst_match=pl.col("match").max().over("name"),
        worst_match_max=pl.col("match_max").max().over("name"),
    )
    return y


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
    check_align(pipe(sam, parse_sam, parse_cigar), allowed, 17)
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
    .sort("bad_tm", "priority")
)
checked


def final_filter(df: pl.DataFrame):
    sortd = df.sort(
        [
            "bad_tm",
            "priority",
            "n_mapped",
        ],
        descending=False,
    ).with_columns(
        final=(pl.col("constructed") + "TATAGTGAGTCGTATTA" + footer),
    )  # ps_ir29
    # sortd["constructed"] += "TATAGTGAGTCGTATTAGACCGGTCT"  # ps_ir49
    return sortd


finale = checked[:66].groupby("gene").apply(final_filter)
# finale = finale.astype(dict(code1=int, code2=int, code3=int, code4=int))
# finale.to_csv("tricycle_probes.csv")

# %%
t7 = reverse_complement("TAATACGACTCACTATAGGG")

assert all(finale["final"].str.contains(t7))

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
# %%# %%
# %%# %%
# %%# %%
# pl.read_csv("all_probes.csv", has_header=False).with_columns(p2=pl.col( "column_2").str.slice(43,20).apply(reverse_complement)).join(zero, left_on="p2", right_on="seq", how="left")['name'].value_counts().sort("name")
