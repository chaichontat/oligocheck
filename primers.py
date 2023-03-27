#%%
import random
from binascii import hexlify
from collections import Counter
from itertools import product

import colorama
import numpy as np
import pandas as pd
import primer3
import seaborn as sns
from hexhamming import hamming_distance_string
from Levenshtein import distance
from nupack import Model, SetSpec, Strand, Tube, tube_analysis

df = pd.read_excel("/Users/chaichontat/Downloads/scirnaseq3PlatePrimers.xlsx", sheet_name=1)
sns.set()

p5 = "AATGATACGGCGACCACCGAGATCTACAC"
truseq1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
p7 = "CAAGCAGAAGACGGCATACGAGAT"
nextera2 = "GTCTCGTGGGCTCGG"
smallrna = "GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
sticky = "CTCACTG"


def gen_plate(name: str, seqs: list[str]):
    wells = [f"{row}{col:02d}" for row in "ABCDEFGH" for col in range(1, 13)]
    return pd.DataFrame(
        {
            "Well Position": wells[: len(seqs)],
            "Name": [name + w for w in wells[: len(seqs)]],
            "Sequence": seqs,
        }
    )


def reverse_complement(seq: str):
    return seq[::-1].translate(str.maketrans("ATCG", "TAGC"))


def printc(seq: str):
    for c in seq:
        if c == "A" or c == "a":
            print(colorama.Fore.GREEN + c, end="")
        elif c == "T" or c == "t":
            print(colorama.Fore.RED + c, end="")
        elif c == "C" or c == "c":
            print(colorama.Fore.BLUE + c, end="")
        elif c == "G" or c == "g":
            print(colorama.Fore.YELLOW + c, end="")
        else:
            print(colorama.Fore.WHITE + c, end="")
    print(colorama.Fore.RESET)


def gc_content(seq: str):
    return (seq.count("G") + seq.count("C")) / len(seq)


def hairpin(seq: str, model: Model) -> float:
    hi = Strand(seq.upper(), "hi")
    t1 = Tube(strands={hi: 1e-7}, name="Tube t1", complexes=SetSpec(max_size=2))
    tube_result = tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    return tube_result["(hi)"].mfe_stack


# %%


dic = {0: "A", 1: "T", 2: "G", 3: "C"}


def gen_from_idx(idx: int, n: int = 10) -> str:
    return dic[idx & 3] + gen_from_idx(idx >> 2, n - 1) if n > 0 else ""


def gen_primer(seq: str, type_: str) -> str:
    if type_ == "nextera":
        return p7 + seq.lower() + nextera2
    if type_ == "smallrna":
        # For TM. Full TruSeq is added in the first PCR
        return p7 + seq.lower() + smallrna[:19]
    raise ValueError("Invalid primer type")


#%%

base = ["".join(x) for x in product("ATGC", repeat=10)]

random.seed(42)
random.shuffle(base)
#%%
my_model = Model(material="dna", celsius=63, sodium=0.3, magnesium=0.003)


def gen_set(start: str, max_gen: int = 400):
    conditions = dict(mv_conc=375, dv_conc=3, dntp_conc=0.2, dna_conc=200)

    out = [
        (
            start,
            hairpin(gen_primer(start, "nextera"), my_model),
            hairpin(gen_primer(start, "smallrna"), my_model),
        )
    ]
    for i, x in enumerate(base):
        if len(out) > max_gen:
            break

        if i % 5000 == 0:
            print(i, len(out))

        seq = "".join(x)
        if seq.endswith("CC"):
            continue
        if "GGGG" in seq or "CCCC" in seq or "AAAA" in seq or "TTTT" in seq:
            continue
        if gc_content(seq) < 0.3 or gc_content(seq) > 0.6:
            continue
        for x in out:
            if distance(seq, x[0]) < 4:
                break
        else:
            # if (tm_next := primer3.calc_hairpin_tm(gen_primer(seq, "nextera"), **conditions)) > 55:
            #     continue
            if (mfenext := hairpin(gen_primer(seq, "nextera"), my_model)) < -0.06:
                continue
            if (mfesm := hairpin(gen_primer(seq, "smallrna"), my_model)) < -0.06:
                continue
            out.append((seq, mfenext, mfesm))

    tms = sorted(out, key=lambda x: -x[1])
    return tms


random.seed(10)
out = gen_set("".join(random.choices("ATGC", k=10)), max_gen=800)
print(len(out))
# [gen_set(x) for x in range(10)]

# final = [x[0] for x in tms if x[1] < 58]
# %%


def calculate_base_diversity(seq: list[str]):
    for i in range(len(seq[0])):
        print({k: v / len(seq) for k, v in Counter([x[i] for x in seq]).items()})


calculate_base_diversity([x[0] for x in out[:96]])

# %%
# %%


# out = 0
# for x in "ATCG":
#     out = (out << 4) | (ord(x) & 3)
# %%
def gen_p5(seq: str):
    return p5 + seq + truseq1[:21]


base = ["".join(x) for x in product("ATGC", repeat=10)]

random.seed(50)
random.shuffle(base)
#%%


def gen_set(seed: int, mismatch: int = 4, max_gen: int = 800):
    print(seed)
    random.seed(seed)

    out = [(temp := "".join(random.choices("ATGC", k=10)), hairpin(gen_p5(temp), my_model))]
    for x in base:
        if len(out) > max_gen:
            break
        seq = "".join(x)
        if seq.endswith("CC"):
            continue
        if "GGGG" in seq or "CCCC" in seq or "AAAA" in seq or "TTTT" in seq:
            continue
        if gc_content(seq) < 0.3 or gc_content(seq) > 0.7:
            continue
        for x in out:
            if distance(seq[:-1], x[:-1]) < mismatch:
                break
        else:
            if (g := hairpin(gen_p5(seq), my_model)) < -0.37:
                continue
            out.append((seq, g))

    else:
        print("not enough")

    return out


p5set = gen_set(50)

# %%
mm = []
my_model = Model(material="dna", celsius=63, sodium=0.3, magnesium=0.003)
for i in range(200):
    if i % 100 == 0:
        print(i)

    for j in range(i + 1, 200):
        hi = Strand(gen_p5(p5set[i][0]), "hi")
        hi2 = Strand(gen_primer(out[j][0], "nextera").upper(), "hi2")
        hi3 = Strand(gen_primer(out[j][0], "smallrna").upper(), "hi3")
        t1 = Tube(strands={hi: 1e-7, hi2: 1e-7, hi3: 1e-7}, name="Tube t1", complexes=SetSpec(max_size=2))
        tube_result = tube_analysis(tubes=[t1], compute=["mfe"], model=my_model)

        # Complexes hi+hi2 and hi+hi3
        conc = max(x[1] for x in list(tube_result.tubes[t1].complex_concentrations.items())[2:4])

        mm.append({"p5": i, "p7": j, "conc": conc})
        mm.append({"p5": j, "p7": i, "conc": conc})

# %%
mm = pd.DataFrame(mm)
ne = {}
for i in range(200):
    ne[i] = mm[mm.p5 == i]
    ne[i] = (ne[i].p7.to_list(), ne[i].conc.to_list())

#%%
# maximum bipartite matching
import networkx as nx

G = nx.Graph()
for i in range(96):
    for j, conc in zip(ne[i][0], ne[i][1]):
        G.add_edge(f"p5-{i:02d}", f"p7-{j:02d}", weight=conc)

m = list(
    nx.bipartite.matching.minimum_weight_full_matching(
        G, {n for n in G.nodes() if n.startswith("p5")}
    ).items()
)
pairings = [(int(x[0][3:]), int(x[1][3:])) for x in m if x[0].startswith("p5")]
pairings = sorted(pairings, key=lambda x: x[0])


# %%


# %%


def umi_generator(template: str, seed: int = 42):
    mapping = {
        "A": "A",
        "T": "T",
        "G": "G",
        "C": "C",
        "N": "ATGC",
        "W": "AT",
        "S": "GC",
        "M": "AC",
        "K": "GT",
        "R": "AG",
        "Y": "CT",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
    }
    iters = [mapping[x] for x in template]
    return ["".join(x) for x in product(*iters)]


base = list(umi_generator("".join(["H"] * 12)))
random.shuffle(base)
bridges = pd.read_csv("./ps_bridges.tsv", sep="\t")

conditions = dict(mv_conc=150, dv_conc=3, dntp_conc=0, dna_conc=100)
idxs = []
for i in range(800):
    if i % 100 == 0:
        print(i)
    o = [
        primer3.calc_hairpin_tm(f"{sticky}{bridges.iloc[i,1][:10]}{x}TTTTTTTT", **conditions)
        for x in base[:1000]
    ]
    # sns.kdeplot(o, bw_adjust=0.1)
    idxs.append((len([x for x in o if x == 0]) / len(o), i))

# %%
rt_selected = sorted(idxs, key=lambda x: (-x[0], x[1]))[:96]
rt_selected
# %%
for i in range(96):
    for j in range(i + 1, 96):
        assert distance(bridges.iloc[rt_selected[i][1], 1], bridges.iloc[rt_selected[j][1], 1]) < 5
# %%


# %%

base = ["".join(x) for x in product("ATGC", repeat=11)]

random.seed(42)
random.shuffle(base)


def gen_lig_primer(seq: str, cut: bool = False) -> str:
    if cut:
        seq = seq[:10]
    return reverse_complement(sticky) + seq.lower() + "T" + truseq1 + reverse_complement(seq).lower()


def hamming(a: bytes, b: bytes) -> int:
    return hamming_distance_string(hexlify(a), hexlify(b)) >> 1


def gen_set(seed: int, mismatch: int = 4, max_gen: int = 800):
    print(seed)
    random.seed(seed)

    out = ["".join(random.choices("ATGC", k=11))]
    for x in base:
        if len(out) > max_gen:
            break
        seq = "".join(x)
        if seq.endswith("CC"):
            continue
        if "GGGG" in seq or "CCCC" in seq or "AAAA" in seq or "TTTT" in seq:
            continue
        if gc_content(seq) < 0.3 or gc_content(seq) > 0.7:
            continue
        for x in out:
            if distance(seq[:-1], x[:-1]) < mismatch:
                break
        else:
            out.append(seq)
    else:
        print("not enough")

    return out


ligs = gen_set(40)
#%%
ligplate = gen_plate("sciv2-Lig-", [gen_lig_primer(x, i & 1) for i, x in enumerate(ligs[:96])])
htvn = "".join(["H"] * 12) + "".join(["T"] * 23) + "VN"
rtplate = gen_plate(
    "sciv2-RT-", ["/5Phos/" + sticky + bridges.iloc[s[1], 1][:10].lower() + htvn for s in rt_selected]
)
finalp5 = [gen_p5(p5set[x[0]][0].lower()) for x in pairings]
finalp7next = [gen_primer(out[x[1]][0], "nextera") for x in pairings]
finalp7small = [gen_primer(out[x[1]][0], "smallrna") for x in pairings]
with pd.ExcelWriter("v2/all.xlsx") as excel:
    gen_plate("sci-v2-P5-", finalp5).to_excel(excel, sheet_name="sci-v2-P5", index=False)
    gen_plate("sci-v2-P7-Nextera-", finalp7next[:24]).to_excel(
        excel, sheet_name="sci-v2-P7-Nextera", index=False
    )
    gen_plate("sci-v2-P7-smallrna-", finalp7small[:24]).to_excel(
        excel, sheet_name="sci-v2-P7-smallrna", index=False
    )
    rtplate.to_excel(excel, sheet_name="sci-v2-RT", index=False)
    ligplate.to_excel(excel, sheet_name="sci-v2-Lig", index=False)
# %%

bridges
# %%
from collections import defaultdict

kmers = defaultdict(int)
for s in bridges.seq:
    for i in range(len(s) - 20):
        kmers[(s[i : i + 20])] += 1
# %%
# 5' AGACGTAAGCACCCTAACGC
# 3' TGTGACTGCTCCTAATACGACTCAC
cutted = bridges.seq.map(lambda x: x[:20])
cutted = pd.concat([cutted, bridges.seq.map(lambda x: x[-20:])]).reset_index(drop=True)
skipped = []
for i in range(len(cutted)):
    for j in range(i + 1, len(cutted)):
        if j in skipped:
            continue
        if distance(cutted[i], cutted[j]) < 7:
            skipped.append(j)
            break
# %%
cutted[cutted.index.difference(skipped)]

#%%
vvv = pd.read_csv("readouts.txt", sep=",", header=None)
vvv[1] = vvv[1].map(lambda x: x.replace(" ", ""))
# %%
vvv["scale"] = "25nm"
vvv["pur"] = "STD"
# %%
vvv.to_csv("readouts2.csv", index=False, header=None, sep="\t")
# %%
