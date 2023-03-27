#%%
from itertools import combinations, cycle, permutations
from typing import Any

import matplotlib.pyplot as plt
import mygene
import numpy as np
import pandas as pd
import primer3
import seaborn as sns
from Bio import SeqIO
from Levenshtein import distance

df = pd.read_csv("/Users/chaichontat/out.csv")


sns.set()


def gen_top(n: int, start: int = 50):
    lenout = 0
    while lenout < n:
        dfs = (
            df.iloc[df["pc1.rot"].map(abs).sort_values(ascending=False).head(start).index],
            df.iloc[df["pc2.rot"].map(abs).sort_values(ascending=False).head(start).index],
        )
        lenout = len(pd.concat(dfs, axis=0)["symbol"].unique())
        start += 5
    return dfs, dfs[0].index.union(dfs[1].index)


dfs, x = gen_top(120)
# %%
res = mygene.MyGeneInfo().getgenes(
    df.iloc[x]["ensembl"].tolist(), fields="symbol,ensembl.transcript", species="mouse"
)
# %%
fpkm = pd.read_csv("/Users/chaichontat/Downloads/fpkm.tsv", sep="\t")
fpkm["transcript_id"] = fpkm["transcript_id"].map(lambda x: x.split(".")[0])
fpkm.set_index("transcript_id", inplace=True)
# %%
def get_fpkm(res: dict[str, Any]):
    out = {}
    total = 0
    notfound = 0
    for r in res:
        ts = r["ensembl"]["transcript"]
        out[r["symbol"]] = {}
        if isinstance(ts, str):
            ts = [ts]
        total += len(ts)
        for t in ts:
            if t in fpkm.index:
                out[r["symbol"]][t] = fpkm.loc[t, "FPKM"]
            else:
                notfound += 1
    assert all(map(len, out.values())), "Not all genes have transcripts"
    return out


out = get_fpkm(res)
# %%
flattened = {t: f for v in out.values() for t, f in v.items()}
flattensum = {gene: sum(v.values()) for gene, v in out.items()}
# %%
plt.hist([np.log10(x + 1) for x in flattensum.values()], bins=100)
# %%
normal_exp = [
    "Ace2",
    "Adora2a",
    "Aldh1l1",
    "Amigo2",
    "Ano3",
    "Aqp4",
    "Ar",
    "Arhgap36",
    "Avpr1a",
    "Avpr2",
    "Baiap2",
    "Bdnf",
    "Bmp7",
    "Brs3",
    "Calcr",
    "Cbln1",
    "Cbln2",
    "Cckar",
    "Cckbr",
    "Ccnd2",
    "Cd24a",
    "Cdkn1a",
    "Cenpe",
    "Chat",
    "Coch",
    "Col25a1",
    "Cplx3",
    "Cpne5",
    "Creb3l1",
    "Crhbp",
    "Crhr1",
    "Crhr2",
    "Cspg5",
    "Cxcl14",
    "Cyp19a1",
    "Cyp26a1",
    "Dgkk",
    "Ebf3",
    "Egr2",
    "Ermn",
    "Esr1",
    "Etv1",
    "Fbxw13",
    "Fezf1",
    "Fn1",
    "Fst",
    "Gabra1",
    "Gabrg1",
    "Gad1",
    "Galr1",
    "Galr2",
    "Gbx2",
    "Gda",
    "Gem",
    "Gjc3",
    "Glra3",
    "Gpr165",
    "Greb1",
    "Grpr",
    "Htr2c",
    "Igf1r",
    "Igf2r",
    "Irs4",
    "Isl1",
    "Kiss1r",
    "Klf4",
    "Lepr",
    "Lmod1",
    "Lpar1",
    "Man1a",
    "Mc4r",
    "Mki67",
    "Mlc1",
    "Myh11",
    "Ndnf",
    "Ndrg1",
    "Necab1",
    "Nos1",
    "Npas1",
    "Npy1r",
    "Npy2r",
    "Ntng1",
    "Ntsr1",
    "Nup62cl",
    "Omp",
    "Onecut2",
    "Opalin",
    "Oprd1",
    "Oprk1",
    "Oprl1",
    "Oxtr",
    "Pak3",
    "Pcdh11x",
    "Pdgfra",
    "Pgr",
    "Plin3",
    "Pnoc",
    "Pou3f2",
    "Prlr",
    "Ramp3",
    "Rgs2",
    "Rgs5",
    "Rnd3",
    "Rxfp1",
    "Scgn",
    "Selplg",
    "Sema3c",
    "Sema4d",
    "Serpinb1b",
    "Serpine1",
    "Sgk1",
    "Slc15a3",
    "Slc17a6",
    "Slc17a7",
    "Slc17a8",
    "Slc18a2",
    "Slco1a4",
    "Sox4",
    "Sox6",
    "Sox8",
    "Sp9",
    "Synpr",
    "Syt2",
    "Syt4",
    "Sytl4",
    "Tacr1",
    "Tacr3",
    "Tiparp",
    "Tmem108",
    "Traf4",
    "Trhr",
    "Ttn",
    "Ttyh2",
]


high_exp = [
    "Oxt",
    "Penk",
    "Sst",
    "Tac1",
    "Gal",
    "Cartpt",
    "Vgf",
    "Trh",
    "Nts",
    "Scg2",
    "Gnrh1",
    "Tac2",
    "Cck",
    "Crh",
    "Ucn3",
    "Adcyap1",
    "Nnat",
    "Sln",
    "Mbp",
    "Th",
    "Fos",
]
norm = mygene.MyGeneInfo().querymany(
    normal_exp, fields="symbol,ensembl.transcript", scopes="symbol", species="mouse"
)
hi = mygene.MyGeneInfo().querymany(
    high_exp, fields="symbol,ensembl.transcript", scopes="symbol", species="mouse"
)
#%%


hisum = {gene: sum(v.values()) for gene, v in get_fpkm(hi).items()}
# %%
normsum = {gene: sum(v.values()) for gene, v in get_fpkm(norm).items()}
# %%
plt.hist([np.log10(x + 1) for x in hisum.values()], bins=100)
plt.hist([np.log10(x + 1) for x in normsum.values()], bins=100)
# %%
import anndata as ad

adata = ad.read("/Users/chaichontat/Downloads/local.h5ad")

coords = adata.obs[adata.obs["donor_id"] == "MsBrainAgingSpatialDonor_1"][["center_x", "center_y"]].rename(
    columns=dict(center_x="x", center_y="y")
)
coords.index.name = "id"
coords.to_csv("coords.csv")
# %%
data = (
    adata[adata.obs["donor_id"] == "MsBrainAgingSpatialDonor_1"]
    .to_df()
    .rename(columns=adata.var["feature_name"].to_dict())
)
data.index.name = "id"

# %%
linear = 2**data
linear /= linear.sum(axis=0)
entropy = np.sum(linear * np.log(linear), axis=0)


# %%
df = pd.read_csv("/Users/chaichontat/Downloads/celltype.csv")
# %%


primer3.calc_tm("CCTGGCTGACAGCTAATACGACTC")
primer3.calc_hairpin("CCTGGCTGACAGCTAATACGACTC")
primer3.calc_homodimer("CCTGGCTGACAGCTAATACGACTC")
# %%

# def get_union(df, genes):
#     return df[df.Gene.isin(genes)].groupby("Gene").sum().sum(axis=1).sort_values(ascending=False)

# %%
readouts = df["Probe Sequence"].map(lambda x: x[19:39]).unique()

readouts = [
    "ATCCTCCTTCAATACATCCC",
    "ACACTACCACCATTTCCTAT",
    "ACTCCACTACTACTCACTCT",
    "ACCCTCTAACTTCCATCACA",
    "ACCACAACCCATTCCTTTCA",
    "TTTCTACCACTAATCAACCC",
    "ACCCTTTACAAACACACCCT",
    "TCCTATTCTCAACCTAACCT",
    "TATCCTTCAATCCCTCCACA",
    "ACATTACACCTCATTCTCCC",
    "TTTACTCCCTACACCTCCAA",
    "TTCTCCCTCTATCAACTCTA",
    "ACCCTTACTACTACATCATC",
    "TCCTAACAACCAACTACTCC",
    "TCTATCATTACCCTCCTCCT",
    "TATTCACCTTACAAACCCTC",
    "AAACACACACTAAACCACCC",
    "AACTCATCTCAATCCTCCCA",
    "TATCTCATCAATCCCACACT",
    "TCTATCATCTCCAAACCACA",
    "TCCAACTCATCTCTAATCTC",
    "AATACTCTCCCACCTCAACT",
    "ATAAATCATTCCCACTACCC",
    "ACCCAACACTCATAACATCC",
    "TTCCTAACAAATCACATCCC",
    "TTCTTCCCTCAATCTTCATC",
    "TACTACAAACCCATAATCCC",
    "TCCTCATCTTACTCCCTCTA",
    "AATCTCACCTTCCACTTCAC",
    "TTACCTCTAACCCTCCATTC",
    "ACTTTCCACATACTATCCCA",
    "ACCTTTCTCCATACCCAACT",
    "TCAAACTTTCCAACCACCTC",
    "ACACCATTTATCCACTCCTC",
    "TCCCAACTAACCTAACATTC",
    "ACATCCTAACTACAACCTTC",
    "ATCCTCACTACATCATCCAC",
    "TCTCACACCACTTTCCTCAT",
    "TCCCTATCAATCTCCATAAC",
    "AACTACATACTCCCTACCTC",
]


def reverse_complement(seq):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement.get(base, base) for base in reversed(seq))


readouts = list(map(reverse_complement, readouts))
# %%


# %%

footer = "TATTTCCCTATAGTGAGTCGTATTAAGATCGGAAGAGCGTCGT"

# %%


def find_gene_name(x: str):
    idx = 3 if len(x.split(" ")[2]) == 6 else 2
    return x.split(" ")[idx].split("__")[0]


# %%
# Use Pclaf Cenpf Cenpa Hells
# %%
ps = pd.read_csv("/Users/chaichontat/Downloads/mm39_refseq_newBalance.tsv", sep="\t", header=None)
#%%
picked2 = ps[ps[13].isin(["Pclaf", "Hells", "Cux2", "Cpne7", "Rapgefl1"])]
picked = df[
    df.Gene.isin(
        [
            "Slc17a7",
            "Gad2",
        ]
    )
]
picked.loc[:, "Probe Sequence"] = picked["Probe Sequence"].apply(lambda x: x[19:112])
picked.rename(columns={"Gene": "gene", "Probe Sequence": "seq"}, inplace=True)
# %%

# %%
# %%


# %%
for i in range(5):
    print(pps.iloc[tmss[i][0], 1])
    primer3.calc_homodimer(pps.iloc[tmss[2][0], 1], **conditions)
    primer3.calc_heterodimer_tm(pps.iloc[tmss[2][0], 1], "CGTATTAAGATCGGAAGAGCGTCGT", **conditions)

# %%
primer3.calc_heterodimer_tm(pps.iloc[tmss[2][0], 1], "CGTATTAAGATCGGAAGAGCGTCGT", **conditions)

tms = []
for seq in tosearch:
    tms.append((seq, primer3.calc_hairpin_tm(pps.iloc[42, 1] + seq, **conditions)))
# %%
primer3.calc_homodimer_tm("TTAATTCAGGGTCCGCGTAGCGAGGGTTTGTAAGGTGAATAATCCTTATGT")

# %%
ps[ps[13] == "Pclaf"][3].apply(len)
# %%
# %%
bridges = pd.read_csv("./ps_bridges.tsv", sep="\t")

# %%


# %%
def slide(x: str, n: int = 20):
    return [x[i : i + n] for i in range(len(x) - n + 1)]


def pick_best(seq: str, target: float = 51):
    return min(slide(seq), key=lambda x: (primer3.calc_tm(x[1], **conditions) - target) ** 2)


dTmarker = "TTACACTCCATCCACTCAA"

for i in range(4):
    print(pick_best(bridges.seq[i]))
# %%


def get_probe_seq(seq: str):
    return seq[19:39], seq[40:70], seq[71:91], seq[92:112]


df.iloc[:100, 1].apply(get_probe_seq)

mhd4_7 = np.array(
    [
        [1, 0, 1, 1, 1, 0, 0],
        [1, 1, 1, 0, 0, 1, 0],
        [0, 1, 0, 1, 1, 1, 0],
        [0, 1, 1, 1, 0, 0, 1],
        [1, 1, 0, 0, 1, 0, 1],
        [1, 0, 0, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 1],
    ]
)
codes = np.apply_along_axis(np.nonzero, axis=1, arr=mhd4_7).squeeze()

#%%

code_dict = {i: pick_best(bridges.seq[i]) for i in range(11)}
code_dict_rev = {i: reverse_complement(pick_best(bridges.seq[i])) for i in range(11)}

mygenes = ["Rapgefl1", "Top2a", "Cux2", "Cpne7"]
selected = ps[ps[13].isin(mygenes)][[13, 3]].rename(columns={13: "gene", 3: "seq"}).reset_index(drop=True)

selected.loc[(selected.gene == "Cux2") & (selected.index % 2), "gene"] = "Cux2_1"
selected.loc[(selected.gene == "Cux2") & ~(selected.index % 2), "gene"] = "Cux2_2"

#%%

out = []
for code, gene in zip(codes, selected.gene.unique()):
    count = 0

    for idxs, seq in zip(cycle(permutations(code, 3)), selected.loc[selected.gene == gene, "seq"]):
        if "CCCC" in seq or "GGGG" in seq:
            continue

        out.append(
            dict(
                gene=gene,
                seq=f"{code_dict_rev[idxs[0]]}A{seq[:30]}A{code_dict_rev[idxs[1]]}A{code_dict_rev[idxs[2]]}",
            )
        )
        count += 1
out = pd.concat([picked, pd.DataFrame(out)])
# %%
conditions = dict(
    mv_conc=300,
    dv_conc=0,
    dntp_conc=0,
    dna_conc=5,
)
# tosearch = combi.seq.map(lambda x: x[:30]).unique()
paintshopprimers = ["ps_ir.tsv", "ps_if.tsv", "ps_or.tsv", "ps_of.tsv"]
pps = pd.concat([pd.read_csv(f"/Users/chaichontat/Downloads/{x}", sep="\t") for x in paintshopprimers])
tmss = []
for i, p in enumerate(pps.seq):
    tms = []
    for seq in out.seq:
        tms.append(primer3.calc_hairpin_tm(p + seq[:30], **conditions))
    tmss.append((i, max(tms)))
    # for seq in out.seq:
    #     tms.append(primer3.calc_hairpin_tm(reverse_complement(p) + seq[:30], **conditions))
    # tmss.append((f"{i}r", max(tms)))

tmss = sorted(tmss, key=lambda x: x[1])


# %%
header = "AGACGTAAGCACCCTAACGC"
res = "TATTTCCCTATAGTGAGTCGTATTAG" + "GAGCAGTCACA"
# %%
outout = out.copy()
outout["seq"] = outout.seq.map(lambda x: header + x + res)
# %%
# outout.to_csv("outout.csv", index=False, header=None)
# %%
q5 = dict(mv_conc=375, dv_conc=2, dntp_conc=0.2, dna_conc=500)
# %%
primer_f = "AGACGTAAGCACCCTAACGC"
primer_r = "TGTGACTGCTCCTAATACGACTCAC"

# %%
# %%
