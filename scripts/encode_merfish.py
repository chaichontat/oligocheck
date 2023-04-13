# %%

import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable

import numpy as np
import numpy.typing as npt
import pandas as pd

from oligocheck.sequtils import equal_distance, gc_content, reverse_complement

wants = [
    "Pclaf",
    "Cenpa",
    "Hells",
    "Top2a",
    "Mcm3",
    "Tipin",
    "Cenpf",
    "Mcm6",
    "Cdc20",
    "Ube2c",
    "Cdca8",
    "Birc5",
    "Smc2",
    "Rpa2",
    "Pcna",
    "Ccne1",
    "Tpx2",
    "Cks2",
    "Mki67",
    "Cdk1",
    "Prc1",
    "Lig1",
    "Ccnb1",
    "Mcm5",
    "Chaf1b",
    "Smc4",
    "Cenpe",
    "Dtl",
    "Uhrf1",
    "Spc25",
    "Nusap1",
    "Racgap1",
    "Chaf1a",
    "Mcm2",
    "Lmnb1",
    "Pbk",
    "Fen1",
    "Knstrn",
    "Mcm7",
    "Gmnn",
    "E2f1",
    "Tacc3",
    "Nasp",
    "Cdkn3",
    "Ccne2",
    "Clspn",
    "Ccnb2",
    "Incenp",
    "Ccng1",
    "Ezh2",
    "Topbp1",
    "Aurka",
    "Cenpq",
    "Kif22",
    "Aurkb",
    "Mcm4",
    "Mis18bp1",
    "Knl1",
    "Stmn1",
    "Ckap2",
    "Dbf4",
    "Psrc1",
    "Rif1",
    "Sgo2a",
    "Cdc6",
    "Fbxo5",
    "Rad21",
    "Ccnd3",
    "Kif4",
    "Sgo1",
    "Timeless",
    "Aspm",
    "Aunip",
    "Exo1",
    "Cenpw",
    "Ube2s",
    "Mad2l1",
    "Hjurp",
    "Kif11",
    "Nuf2",
    "Pmf1",
    "Cdkn2c",
    "Kif20b",
    "Esco2",
    "Cdc45",
    "Casp8ap2",
    "Nde1",
    "Cdca2",
    "Atad5",
    "Skp2",
    "Cdc25c",
    "Psmc3ip",
    "Tfdp1",
    "Cdkn1a",
    "Rpa1",
    "Ncapg",
    "Ckap5",
    "Ncapd3",
    "Bub3",
    "Kif2c",
    "Etaa1",
    "Vrk1",
    "Bub1b",
    "Cdk2",
    "Ncapg2",
]


def count_genes(df: pd.DataFrame) -> pd.Series:
    if "gene" not in df.columns:
        df["gene"] = df["name"].apply(lambda x: x.split("_")[0])
    return df.groupby("gene")["gene"].count().sort_values(ascending=False)


def count_oks(df: pd.DataFrame):
    oks = df.columns[df.columns.str.contains("ok_")]
    df["oks"] = df[oks].sum(axis=1)
    return df[df.oks > 3]


# %%
res = []
for gene in wants:
    # try:
    res.append(pd.read_parquet(f"output/{gene}_high.parquet"))
    # except FileNotFoundError:
    #     ...
df = pd.concat(res)
df["gene"] = df["name"].apply(lambda x: x.split("_")[0])
df_ok = count_oks(df)


# %%
def cSpecStackFilter(seq: str):
    seq = reverse_complement(seq)
    for i in range(6):
        if seq[i : i + 6].count("C") >= 4:
            return False
    return True


def ACompfilter(seq: str):
    return seq.count("T") / len(seq) < 0.28


def mark_ok(df: pd.DataFrame):
    # https://www.nature.com/articles/s41596-022-00750-2
    df["ok_quad_c"] = ~df["seq"].str.contains("GGGG")
    df["ok_quad_a"] = ~df["seq"].str.contains("TTTT")
    df["ok_comp_a"] = df["seq"].map(ACompfilter)
    df["ok_stack_c"] = df["seq"].map(cSpecStackFilter)
    df["gc_content"] = df["seq"].map(gc_content)
    df["ok_gc"] = (df.gc_content >= 0.35) & (df.gc_content <= 0.65)
    df["ok_hairpin"] = df.hp < 34
    return df


def handle_done(df: pd.DataFrame, target: int = 40):
    if len(df) < target:
        return df
    picked = df.sort_values(by="pos", ascending=False)
    return picked.iloc[equal_distance(len(picked), target)]


def handle_overlap(tss: Iterable[str], filtered: pd.DataFrame, overlap: int = -1):
    selected = set()
    for ts in tss:
        nool = sorted(
            [(r, idx, pos) for idx, r in filtered.iterrows() if (pos := get_seq(ts).find(r.seq)) > -1],
            key=lambda x: x[2],
        )
        curr = min(0, -overlap)
        if len(nool[1][0].seq) < overlap:
            raise ValueError("Overlap too large")
        for r, idx, pos in nool:
            if pos - curr >= 0:
                selected.add(idx)
                curr = pos + len(r.seq) - overlap
    return filtered.loc[list(selected)]


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


# pick(df[df.gene == "Aspm"])
df = mark_ok(df)
filtered = count_genes(df[df.oks > 4])
len(filtered[filtered > 40].index)

pick(df[df.gene == "Aspm"]).sort_values(by="pos")

# %%
lessthan50 = counts[counts < 50].index
levels = ["high", "medium", "low", "high_ol", "medium_ol", "low_ol"]


def runpls(gene: str):
    # if not Path(f"output/{gene}_low.parquet").exists():
    print("running", gene)
    subprocess.run(["python", "samwhat.py", gene, "low_ol", "--output", "output/"])


with ThreadPoolExecutor(12) as executor:
    executor.map(runpls, lessthan50)
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
