# %%


import re
import subprocess
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from itertools import cycle, permutations
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
import numpy.typing as npt
import pandas as pd
import polars as pl
import primer3

from oligocheck.merfish.encoding import filter_specifity
from oligocheck.merfish.external_data import all_transcripts, gene_to_eid, get_gencode, get_seq
from oligocheck.merfish.filtration import the_filter
from oligocheck.sequtils import (
    equal_distance,
    formamide_correction,
    gc_content,
    parse_sam,
    plot_local_gc_content,
    reverse_complement,
    tm_hybrid,
)

wants = pd.read_csv("tricyclegenes.csv", header=None)[0][:110]


def count_genes(df: pd.DataFrame) -> pd.Series:
    if "gene" not in df.columns:
        df["gene"] = df["name"].apply(lambda x: x.split("_")[0])
    return df.groupby("gene")["gene"].count().sort_values(ascending=False)


def filter_gene(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    return df[df["gene"] == gene]


def handle_done(df: pd.DataFrame, target: int = 40):
    if len(df) < target:
        return df
    picked = df.sort_values(by="pos", ascending=False)
    return picked.iloc[equal_distance(len(picked), target)]


everything = get_gencode("data/mm39/gencode_vM32_transcripts.parquet")


def count_oks(df: pd.DataFrame, max_notok: int = 2):
    total_oks = len(df.columns[df.columns.str.contains("ok_")])
    oks = df.columns[df.columns.str.contains("ok_")]
    df["oks"] = df[oks].sum(axis=1)
    return df[df.oks > total_oks - max_notok]


# %%
temp = []
for gene in wants:
    try:
        temp.append(pd.read_parquet(f"output/{gene}_medium.parquet"))
    except FileNotFoundError:
        print("File not found", gene)

df = pd.concat(temp)
total_oks = len(df.columns[df.columns.str.contains("ok_")])
oks = df.columns[df.columns.str.contains("ok_")]
df["oks"] = df[oks].sum(axis=1)
# df = mark_ok(df)


# %%


# %%
with ThreadPoolExecutor(8) as executor:
    out = list(executor.map(lambda x: the_filter(df, x), np.array_split(wants, 8)))
out = pd.concat(out).set_index("name")
# %%

# eid = gene_to_eid(gene)
# tss_gencode = tuple(everything[everything.gene_id == eid]["transcript_id"].values)
# tss_all = all_transcripts(gene)
# out = pd.concat(out)
counts = count_genes(out)
easy = counts[counts >= 50].index
noteasy = counts[counts < 50].index
res = the_filter(df, noteasy, overlap=15).set_index("name")

out = pd.concat([out[out.gene.isin(easy)], res])
# %%


# out = pd.concat([out, res])

# %%
wants_filtered = [x for x in wants if (x != "Stmn1" and x not in noteasy[-9:])]
# %%

# %%
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 5), sharex=True)
sns.stripplot(out, x="pos_start", alpha=0.5, ax=ax1)
plot_local_gc_content(ax2, get_seq(filter_gene(df, "Timeless").transcript_ori.unique()[0]), 50, alpha=1)


# %%
fpkm = pd.read_csv("fpkm.tsv", sep="\t", index_col=0)
fpkm.index = fpkm.index.str.split(".").str.get(0)

fp = {}
for gene in list(wants.values) + ["Eef2"]:
    try:
        g = everything[everything.gene_id == gene_to_eid(gene)]
    except KeyError:
        print("not found", gene)
        continue
    g = g[g.transcript_support_level.isin(["1", "2"])]
    tss = g.transcript_id.values
    # tss = all_transcripts(gene)
    s = 0
    for ts in tss:
        try:
            s += fpkm.loc[ts, "FPKM"]
        except KeyError:
            print(gene, ts)
    fp[gene] = s
sorted(fp.items(), key=lambda x: x[1])
# %%
# %%
fpfil = [fp[x] for x in wants_filtered]


# %%
mhd4 = np.loadtxt("./140gene16bit.csv", delimiter=",", dtype=bool)


def gen_codebook(seed: int):
    rand = np.random.RandomState(seed)
    rmhd4 = mhd4.copy()
    rand.shuffle(rmhd4)
    return rmhd4


def find_optimalish(seed: int, fpkm: Iterable[float]):
    rmhd4 = gen_codebook(seed)
    res = rmhd4[:100] * np.array(fpkm)[:, np.newaxis]
    tocalc = res.sum(axis=0)
    normed = tocalc / tocalc.sum()
    return -np.sum(normed * np.log2(normed))


k = [find_optimalish(i, fpfil) for i in range(1000)]
best = np.argmax(k)

y = gen_codebook(best)
(y[:100] * np.array(fpfil)[:, np.newaxis]).sum(axis=0)
# %%
# %%
zeroth_readout = pd.read_csv("zeroth_readouts.csv", index_col=0)
fusedreadout = []
for i in zeroth_readout.index:
    for j in zeroth_readout.index:
        if i == j:
            continue
        fusedreadout.append(
            dict(name=i + j, seq="AA" + zeroth_readout.loc[i] + "AA" + zeroth_readout.loc[j] + "AA")
        )


fusedreadout = pd.DataFrame(fusedreadout).set_index("name")["seq"]

with open("zeroth.fastq", "w") as f:
    for i, row in zip(fusedreadout.index, fusedreadout):
        f.write(f"@{i}\n")
        f.write(row.seq + "\n")
        f.write("+\n")
        f.write("~" * len(row.seq) + "\n")
# %%
zero = parse_sam("zeroth.sam", split_name=False)
zero = check_align(zero)
zero["bad_tm"] = zero.tm_offtarget > 33

blacklist = set(
    map(
        tuple,
        (
            zero[zero.bad_tm]
            .name.str.split(r"Readout-(\d+)-v2Readout-(\d+)-v2", expand=True)
            .iloc[:, 1:3]
            .astype(int)
            .to_numpy()
            .tolist()
        ),
    )
)


# %%


def gen_full_probe(
    df: pd.DataFrame,
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
            raise ValueError("Not enough probes for gene", gene)
        code_ids = np.argwhere(codeset).squeeze()
        codeseqs = zeroth_readout.iloc[code_ids].to_numpy().squeeze()
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
res = "TATTTCCCTATAGTGAGTCGTATTAG" + "GGACACCTACAG"

for_bowtie = gen_full_probe(out, easy, y, header, "TATTTCC")
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

bridges = pd.read_csv("ps_bridges.tsv", sep="\t")
with open("bridges.fastq", "w") as f:
    for i, row in bridges.iterrows():
        f.write(f"@{row['id']}\n")
        f.write(row.seq + "\n")
        f.write("+\n")
        f.write("~" * len(row.seq) + "\n")


# %%
sam = parse_sam("final.sam")

# %%

c = re.compile(r"(\d+)S(\d+)M")
c2 = re.compile(r"(\d+)M")


def parse_cigar(s: str):
    try:
        return list(map(int, c.findall(s)[0]))
    except IndexError:
        try:
            return 0, int(c2.findall(s)[0])
        except IndexError:
            return 0, 0


# %%


# %%
# out = []
# for i, rows in sam.groupby("gene"):
#     out.append(filter_specifity(all_transcripts(rows.iloc[0]["gene"]), rows))

# %%


sam["align_start"], sam["matches"] = zip(*sam.cigar.map(parse_cigar))
# set dtype of matches to int
# sam = sam.astype(dict(matches=int, align_start=int))
to_check = sam[sam["matches"] > 17]


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


offtarget = pd.concat(check_align(to_check, gene).sort_values("tm_offtarget") for gene in easy).set_index(
    "name"
)

# %%
# for_bowtie.set_index("name", inplace=True)
joined = for_bowtie.join(offtarget[["tm_offtarget"]], how="outer").fillna(0)
joined["bad_tm"] = joined.tm_offtarget > 35
# for rows in joined.groupby("gene"):
joined[joined.gene == "Mki67"].sort_values(
    ["bad_tm", "n_mapped", "priority"], ascending=[True, True, False]
).iloc[:40]
# %%
# %%
