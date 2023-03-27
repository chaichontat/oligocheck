#%%
from itertools import chain

import pandas as pd
import primer3
from Levenshtein import distance

from oligocheck.utils import reverse_complement

conditions = dict(
    mv_conc=300,
    dv_conc=0,
    dntp_conc=0,
    dna_conc=3,
)
dTmarker = "TTACACTCCATCCACTCAA"


def slide(x: str, n: int = 20) -> list[str]:
    return [x[i : i + n] for i in range(len(x) - n + 1)]


def pick_best(seq: str, target: float = 51):
    return (
        min(slide(seq)[:24], key=lambda x: (primer3.calc_tm(x[1], **conditions) - target) ** 2),
        min(slide(seq)[22:], key=lambda x: (primer3.calc_tm(x[1], **conditions) - target) ** 2),
    )


bridges = pd.read_csv("./ps_bridges.tsv", sep="\t")
cutted = list(chain.from_iterable([pick_best(x.seq) for x in bridges.itertuples()]))

#%%
selected = [dTmarker]
generated = []
dist = 6
for seq in bridges.itertuples():
    for cutted in slide(seq.seq):
        if distance(cutted, selected[-1]) < dist:
            continue
        if not (49 <= primer3.calc_tm(cutted, **conditions) <= 53):
            continue

        for s in selected:
            if distance(s, cutted) < dist:
                break
        else:
            if (
                primer3.calc_hairpin_tm(
                    gen := "AA" + reverse_complement(cutted) + "TAAAAAAAAAAAAAAAAAAAAAA", **conditions
                )
                <= 0
            ):
                selected.append(cutted)
                generated.append(gen)
# %% Prep for IDT
def genidt(seq: str) -> str:
    return "/5AmMC6/" + seq[:-2] + "*A*A"


idted = [genidt(x) for x in generated]
#%%

proteins = ["lectin", "mouse", "rabbit", "rat"]
for i, (name, sel, conj) in enumerate(zip(proteins, reversed(selected), reversed(idted)), 1):
    print(f"MERConjRev-{i:04d}-{name}\t{conj}\t\t100nm\tSTD")
    print(f"ReadoutRev-{i:04d}\t/5AmMC6/{sel}\t\t100nm\tSTD")
print()
#%%
