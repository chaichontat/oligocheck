# %%

import subprocess
from concurrent.futures import ThreadPoolExecutor
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
from oligocheck.sequtils import (
    equal_distance,
    gc_content,
    parse_sam,
    plot_local_gc_content,
    reverse_complement,
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


def handle_overlap(
    filtered: pd.DataFrame, criteria: list[Callable[[pd.Series], bool]] = [], overlap: int = -1, n: int = 200
):
    if not filtered.gene.unique().size == 1:
        raise ValueError("More than one gene in filtered")
    filtered = filtered.sort_values(
        by=["is_ori_seq", "transcript_ori", "pos_start", "tm"], ascending=[False, True, True, False]
    )
    eid = gene_to_eid(filtered.gene.iloc[0])
    tss = tuple(everything[everything.gene_id == eid]["transcript_id"].values)
    if not criteria:
        criteria = [lambda _: True]

    selected = []
    for ts in tss:
        this_transcript = filtered[filtered["transcript_ori"] == ts]
        if len(this_transcript) == 0:
            print("No match found for", ts)
            continue
        forbidden = np.zeros(this_transcript.pos_end.max() + 1 + max(0, overlap), dtype=bool)
        for criterion in criteria:
            if len(selected) >= n:
                break
            for idx, r in this_transcript.iterrows():
                if idx in selected or not criterion(r):
                    continue
                if np.any(forbidden[r.pos_start - 1 : r.pos_end + 1 - overlap]):
                    continue
                selected.append(idx)
                forbidden[r.pos_start - 1 : r.pos_end + 1 - overlap] = 1

    return filtered.loc[selected]


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
def the_filter(df: pd.DataFrame, genes: list[str], overlap: int = -1):
    out = []
    for gene in genes:
        out.append(
            handle_overlap(
                df[(df.gene == gene)],
                criteria=[
                    lambda x: (x.tm > 49)
                    & (x.tm < 54)
                    & (x.oks > 4)
                    & (x.hp < 35)
                    & (x.nonspecific_binding < 0.001),
                    lambda x: (x.tm > 49)
                    & (x.tm < 54)
                    & (x.oks > 4)
                    & (x.hp < 40)
                    & (x.nonspecific_binding < 0.01),
                    lambda x: (x.tm > 49)
                    & (x.tm < 54)
                    & (x.oks > 3)
                    & (x.hp < 35)
                    & (x.nonspecific_binding < 0.001),
                    lambda x: (x.tm > 47) & (x.tm < 56) & (x.hp < 40) & (x.oks > 3),
                    lambda x: (x.tm > 46) & (x.tm < 56) & (x.hp < 40) & (x.oks > 3),
                ],
                overlap=overlap,
                n=40,
            )
        )
    return pd.concat(out)


# %%
with ThreadPoolExecutor(8) as executor:
    out = list(executor.map(lambda x: the_filter(df, x), np.array_split(wants, 8)))
out = pd.concat(out)
# %%

# eid = gene_to_eid(gene)
# tss_gencode = tuple(everything[everything.gene_id == eid]["transcript_id"].values)
# tss_all = all_transcripts(gene)
# out = pd.concat(out)
counts = out.groupby("gene")["gene"].count().sort_values(ascending=False)
counts
# %%
easy = counts[counts >= 40]

# %%

noteasy = counts[counts < 40]
res = the_filter(df, noteasy.index, overlap=15)

count_genes(res)
# %%
wants_filtered = [x for x in wants if (x != "Stmn1" and x not in noteasy.index[-9:])]
# %%

# %%
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 5), sharex=True)
sns.stripplot(out, x="pos_start", alpha=0.5, ax=ax1)
plot_local_gc_content(ax2, get_seq(filter_gene(df, "Timeless").transcript_ori.unique()[0]), 50, alpha=1)


# %%
def pick(df: pd.DataFrame, target: int = 40):
    total_oks = len(list(filter(lambda x: "ok_" in x, df.columns)))
    # picked = df[(df.n_mapped == df.n_mapped.max()) & (df.oks == 4)]
    picked = pd.DataFrame()
    for nm in reversed(range(1, df.n_mapped.max() + 1)):
        for oks in reversed(range(total_oks - 2, total_oks + 1)):
            new = df[(df.n_mapped == nm) & (df.oks == oks)]
            print(f"n_mapped {nm} {oks} got {len(new)}")
            picked = pd.concat([picked, new])
            if len(picked) >= target:
                return handle_done(picked, target)
    print("not enough")
    return handle_done(picked, target)


filtered = count_genes(df[df.oks > 4])
len(filtered[filtered > 40].index)

pick(df[df.gene == wants[1]]).sort_values(by="pos")


# %%
df_ol = pd.concat([pd.read_parquet(f"output/{gene}_medium_ol.parquet") for gene in wants])
# %%
df_low = pd.concat([pd.read_parquet(f"output/{gene}_low.parquet") for gene in lessthan50])

# %%
# %%
ncr = pd.read_csv("../mm39/trcombi.txt", sep=" ", header=None, names=["counts"], index_col=0)["counts"]


# %%
def slide(x: str, n: int = 20):
    return [x[i : i + n] for i in range(len(x) - n + 1)]


for seq in df.seq:
    for x in slide(seq, n=14):
        if x in ncr:
            print(x)
            break
# %%
from Bio import SeqIO

f = SeqIO.parse("../mm39/Mus_musculus.GRCm39.ncrna.fa", "fasta")
# %%
for i in zip(f, range(100)):
    print(i)
# %%
f = SeqIO.parse("../mm39/Mus_musculus.GRCm39.ncrna.fa", "fasta")
# %%
for i in zip(f, range(100)):
    print(i)


# %%
def tm(s: str) -> float:
    """Approximately the same as the results from the MATLAB script"""
    return primer3.calc_tm(s, mv_conc=390, dv_conc=0, dna_conc=1) + 4


# %%
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

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


def find_optimalish(seed: int):
    rmhd4 = gen_codebook(seed)
    res = rmhd4[:100] * np.array(fpfil)[:, np.newaxis]
    tocalc = res.sum(axis=0)
    normed = tocalc / tocalc.sum()
    # entropy
    return -np.sum(normed * np.log2(normed))


k = [find_optimalish(i) for i in range(1000)]
best = np.argmax(k)
# %%
y = gen_codebook(best)
(y[:100] * np.array(fpfil)[:, np.newaxis]).sum(axis=0)
# %%
zeroth_readout = pd.read_csv("zeroth_readouts.csv", index_col=0)
# %%

final = []
for gene, codeset in zip(wants_filtered, y):
    probes = out[out.gene == gene]
    if len(probes) < 40:
        continue
    codeseqs = zeroth_readout.iloc[np.argwhere(codeset).squeeze()].values.squeeze()
    constructed = []
    for (i, probe), codes in zip(probes.iterrows(), cycle(permutations(codeseqs, 4))):
        constructed.append(
            (
                # reverse_complement(probe.seq),
                reverse_complement(
                    codes[0] + "T" + codes[1] + "T" + probe.seq + "T" + codes[2] + "T" + codes[3]
                )
            )
        )
    out.loc[out.gene == gene, "constructed"] = constructed

    # final.append(probes.iloc[equal_distance(len(probes), 40)])
out.dropna(inplace=True)
out["constructed"] = out["constructed"].apply(lambda x: header + x + res)

# %%
header = "TGGCGACTGAAGTACGAGTCC"  # ps_or8

res = "TATTTCCCTATAGTGAGTCGTATTAG" + "GGACACCTACAG"
# %%
primer3.calc_heterodimer_tm("GTGAGTCGTATTAG GCCTGTAGAGTGAT", header, **q5)
# %%
# %%
with open("out.fastq", "w") as f:
    for i, row in out.iterrows():
        f.write(f"@{row['name']}\n")
        f.write(row.constructed + "\n")
        f.write("+\n")
        f.write("~" * len(row.constructed) + "\n")

# %%
sam = parse_sam("final.sam")

# %%
import re

c = re.compile(r"(\d+)M")
# %%
sam["matches"] = sam.cigar.map(lambda x: int(c.findall(x)[0]))
# %%
filter_specifity()
