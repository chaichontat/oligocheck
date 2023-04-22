# %%
import subprocess
from concurrent.futures import ThreadPoolExecutor
from functools import reduce
from itertools import cycle, permutations
from typing import Iterable

import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl

from oligocheck.merfish.external_data import all_transcripts, gene_to_eid, get_gencode

# from oligocheck.merfish.filtration import the_filter
from oligocheck.sequtils import equal_distance, parse_cigar, parse_sam, reverse_complement, tm_hybrid

wants = ["Pcna", "Eef2"]
# wants = ["Eef2"]


def count_genes(df: pl.DataFrame) -> pl.DataFrame:
    return df.groupby("gene").count().sort("count")


def filter_gene(df: pl.DataFrame, gene: str) -> pl.DataFrame:
    return df.filter(pl.col("gene") == gene)


# def handle_done(df: p.DataFrame, target: int = 40):
#     if len(df) < target:
#         return df
#     picked = df.sort_values(by="pos", ascending=False)
#     return picked.iloc[equal_distance(len(picked), target)]


everything = get_gencode("data/mm39/gencode_vM32_transcripts.parquet")

# %%
temp = []
for gene in wants:
    try:
        temp.append(pl.read_parquet(f"output/{gene}_medium.parquet"))
    except FileNotFoundError:
        print("File not found", gene)

df = pl.concat(temp)
df = df.with_columns(
    [
        pl.sum(pl.col("^ok_.*$")).alias("oks"),
    ]
)
total_oks = len(list(filter(lambda x: x.startswith("ok_"), df.columns)))


def handle(df: pl.DataFrame, gene: str, criteria: list[pl.Expr], overlap: int = -1, n: int = 200):
    if len(df.select(pl.col("gene").unique())) > 1:
        raise ValueError("More than one gene in filtered")
    df = df.sort(
        by=["is_ori_seq", "transcript_ori", "pos_start", "tm"], descending=[True, False, False, True]
    )
    eid = gene_to_eid(gene)
    tss = tuple(everything[everything.gene_id == eid].index)

    if not criteria:
        criteria = [pl.all("*")]

    df = df.with_columns(
        [
            pl.lit(0, dtype=pl.UInt8).alias("priority"),
            pl.arange(0, len(df), dtype=pl.UInt32).alias("index"),
        ]
    )
    selected = set()

    for ts in tss:
        this_transcript = df.filter(pl.col("transcript_ori") == ts)
        if len(this_transcript) == 0:
            print("No match found for", ts)
            continue

        forbidden = np.zeros(df.select(pl.col("pos_end")).max().item() + 1 + max(0, overlap), dtype=bool)
        priority = 1
        for criterion in criteria:
            if len(selected) >= n:
                break
            for idx, st, end in (
                df.filter(criterion & ~pl.col("index").is_in(selected))
                .select(["index", "pos_start", "pos_end"])
                .iter_rows()
            ):
                if np.any(forbidden[st - 1 : end + 1 - overlap]):
                    continue
                selected.add(idx)
                df[idx, "priority"] = priority
                forbidden[st - 1 : end + 1 - overlap] = 1
            priority += 1

    return df.filter(pl.col("priority") > 0)


def the_filter(df: pl.DataFrame, genes: Iterable[str], overlap: int = -1):
    out = []
    for gene in genes:
        out.append(
            handle(
                df.filter(pl.col("gene") == gene),
                gene,
                criteria=[
                    # fmt: off
                    pl.col("tm").is_between(49, 54) & (pl.col("oks")>4) & (pl.col("hp") < 35) & (pl.col("nonspecific_binding") < 0.001),
                    pl.col("tm").is_between(49, 54) &  (pl.col("oks")>4) & (pl.col("hp") < 40) & (pl.col("nonspecific_binding") < 0.05),
                    pl.col("tm").is_between(47, 56) &  (pl.col("oks")>3) & (pl.col("hp") < 40),
                    pl.col("tm").is_between(46, 56) &  (pl.col("oks")>3) & (pl.col("hp") < 40),
                    pl.col("tm").is_between(46, 56) &  (pl.col("oks")>2) & (pl.col("hp") < 40),
                    # fmt: on
                ],
                overlap=overlap,
                n=100,
            )
        )
    return pl.concat(out)


# %%
out_nool = the_filter(df, wants)
# %%
counts = count_genes(out_nool)
easy = counts.filter(pl.col("count") >= 45)["gene"]
noteasy = counts.filter(pl.col("count") < 45)["gene"]
res = the_filter(df, noteasy, overlap=15)
out = pl.concat([out_nool.filter(pl.col("gene").is_in(easy)), res])
# %%
wants_filtered = wants
# wants_filtered = [x for x in wants.iloc[:110] if (x != "Stmn1") and x not in noteasy[-9:]]
# out = out[out.gene.isin(wants_filtered)]
# %%

# %%
fpkm = pl.read_parquet("data/fpkm/P0_combi.parquet")


fp = {}
for gene in list(wants):
    try:
        g = everything[everything.gene_id == gene_to_eid(gene)]
    except KeyError:
        print("not found", gene)
        continue
    g = g[g.transcript_support_level.isin(["1", "2"])]
    tss = g.index
    # tss = all_transcripts(gene)
    s = 0
    for ts in tss:
        try:
            s += fpkm.filter(pl.col("transcript_id(s)") == ts)["FPKM"].item()
        except KeyError:
            print(gene, ts)
    fp[gene] = s
sorted(fp.items(), key=lambda x: x[1])
# %%
# %%
fpfil = [fp[x] for x in wants_filtered]


# %%
mhd4 = np.loadtxt("./mhd4_16bit.csv", delimiter=",", dtype=bool)


def gen_codebook(seed: int):
    rand = np.random.RandomState(seed)
    rmhd4 = mhd4.copy()
    rand.shuffle(rmhd4)
    return rmhd4


def find_optimalish(seed: int, fpkm: Iterable[float]):
    rmhd4 = gen_codebook(seed)
    res = rmhd4[: len(fpkm)] * np.array(fpkm)[:, np.newaxis]
    tocalc = res.sum(axis=0)
    normed = tocalc / tocalc.sum()
    return -np.sum(normed * np.log2(normed))


k = [find_optimalish(i, fpfil) for i in range(5000)]
best = np.argmax(k)

y = gen_codebook(best)
print(best)
(y[: len(fpfil)] * np.array(fpfil)[:, np.newaxis]).sum(axis=0)
# %%
# %%
zeroth_readout = pd.read_csv("zeroth_readouts.csv", index_col=0)
fusedreadout = []
for i in zeroth_readout.index:
    for j in zeroth_readout.index:
        if i == j:
            continue
        fusedreadout.append(
            dict(name=i + "_" + j, seq="AA" + zeroth_readout.loc[i] + "AA" + zeroth_readout.loc[j] + "AA")
        )


fusedreadout = pd.DataFrame(fusedreadout).set_index("name")["seq"]

with open("zeroth.fastq", "w") as f:
    for i, row in zip(fusedreadout.index, fusedreadout):
        f.write(f"@{i}\n")
        f.write(row.seq + "\n")
        f.write("+\n")
        f.write("~" * len(row.seq) + "\n")

subprocess.run(
    "bowtie2 -x data/mm39/mm39 zeroth.fastq "
    "--no-hd -t -a --local -D 20 -R 3 "
    "-N 0 -L 16 -i C,2 --score-min G,1,6 -S zeroth.sam",
    shell=True,
    check=True,
)
# %%


def calc_tm(row: pd.Series):
    seq = row["seq"][row["align_start"] : row["align_start"] + row["matches"]]
    return tm_hybrid(seq)


def check_align(df: pd.DataFrame, gene: str | None = None):
    if gene is not None:
        tss = all_transcripts(gene)
        check = df[df.gene == gene]
        check = check[~check.transcript.isin(tss)].copy()
    else:
        check = df.copy()

    check["align_start"], check["matches"] = zip(*check.cigar.map(parse_cigar))
    check["tm_offtarget"] = check.apply(calc_tm, axis=1)
    worst_case = (
        check.sort_values("matches").drop_duplicates("name", keep="last").sort_values(["tm_offtarget"])
    )
    # assert worst_case.iloc[40].tm < 40
    return worst_case


zero = parse_sam("zeroth.sam", split_name=False)
zero = zero[zero.transcript != "*"]
zero = check_align(zero)
zero["bad_tm"] = zero.tm_offtarget > 33

blacklist = set(
    map(
        tuple,
        (zero[zero.bad_tm].name.str.split("_", expand=True).to_numpy().tolist()),
    )
)


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

header = "TGGCGACTGAAGTACGAGTCC"  # ps_or8
# header = "GAGAGGCGAGGACACCTACAG"  # ps_or2

# from 18 for cellcycle
for_bowtie = gen_full_probe(out, zeroth_readout, wants_filtered, y, header, "TATTTCCC")
# %%
with open("out.fastq", "w") as f:
    for i, row in for_bowtie.iterrows():
        f.write(f"@{i}\n")
        f.write(row.constructed + "\n")
        f.write("+\n")
        f.write("~" * len(row.constructed) + "\n")

subprocess.run(
    "bowtie2 -x data/mm39/mm39 out.fastq "
    "--no-hd -t -a --local -D 20 -R 3 "
    "-N 0 -L 17 -i C,2 --score-min G,1,6 -S final.sam",
    shell=True,
    check=True,
)


# %%
sam = parse_sam("final.sam")

# %%


sam["align_start"], sam["matches"] = zip(*sam.cigar.map(parse_cigar))
# set dtype of matches to int
# sam = sam.astype(dict(matches=int, align_start=int))
to_check = sam[sam["matches"] > 16]


offtarget = pd.concat(
    check_align(to_check, gene).sort_values("tm_offtarget") for gene in wants_filtered
).set_index("name")

# %%
# for_bowtie.set_index("name", inplace=True)
joined = for_bowtie.join(offtarget[["tm_offtarget"]], how="outer").fillna(0)
joined["bad_tm"] = joined.tm_offtarget > 33
joined["horrible_tm"] = joined.tm_offtarget > 37
finale = []
for i, rows in joined[~joined.horrible_tm].groupby("gene"):
    sortd = rows.sort_values(
        [
            "bad_tm",
            "priority",
            "n_mapped",
        ],
        ascending=True,
    ).iloc[:40]
    sortd["constructed"] += "TATAGTGAGTCGTATTAGAGGCACTG"  # ps_ir29
    # sortd["constructed"] += "TATAGTGAGTCGTATTAGACCGGTCT"  # ps_ir49
    finale.append(sortd)
finale = pd.concat(finale)
finale = finale.astype(dict(code1=int, code2=int, code3=int, code4=int))
finale.to_csv("tricycle_probes.csv")

# %%
t7 = reverse_complement("TAATACGACTCACTATAGGG")

assert all(finale.constructed.str.find(t7) != -1)

codes = set()
for i, row in finale.groupby("gene")[["code1", "code2", "code3", "code4"]].agg(["unique"]).iterrows():
    assert len(row) == 4
    gene_code = tuple(sorted(reduce(lambda x, y: x & y, map(set, row.values))))
    codes.add(gene_code)

assert len(codes) == finale.groupby("gene").count().shape[0]


# %%
count_genes(finale)
# %%
# %%# %%

# %%
count_genes(finale)
# %%
# %%# %%
# %%
# %%# %%

# %%
count_genes(finale)
# %%
# %%# %%
