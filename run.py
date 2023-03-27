#%%
from itertools import combinations
from pathlib import Path

from oligocheck.io import get_whitelist
from oligocheck.tenx import gen_oligo
from oligocheck.utils import hamming

whitelist = get_whitelist()
bcs = list(range(1, 26))
for x, y in combinations(bcs, 2):
    assert hamming(whitelist[x - 1], whitelist[y - 1]) > 4, f"{x} {y} distance less than 5."


def builder(seq: str, i: int) -> str:
    return "/5Biosg/" + seq


def namer(seq: str, i: int) -> str:
    return f"Hash-sci-{i}-nextera"


Path("out.txt").write_text(gen_oligo(bcs, builder=builder, namer=namer, provider="sci"))

# %%
import pandas as pd

wl = get_whitelist()

df = pd.DataFrame(
    [
        {"barcode": wl[0], "name": "SA"},
        {"barcode": wl[1], "name": "CFSE"},
        {"barcode": wl[8], "name": "Hash2"},
        {"barcode": wl[9], "name": "Hash3"},
        {"barcode": wl[10], "name": "Hash4"},
        {"barcode": wl[11], "name": "Hash5"},
        {"barcode": wl[12], "name": "Hash6"},
        {"barcode": wl[13], "name": "Hash7"},
        {"barcode": wl[14], "name": "Hash8"},
        {"barcode": wl[15], "name": "Hash9"},
        {"barcode": wl[16], "name": "Ladder1"},
        {"barcode": wl[17], "name": "Ladder2"},
        {"barcode": wl[18], "name": "Ladder3"},
        {"barcode": wl[19], "name": "Ladder4"},
        {"barcode": wl[20], "name": "Ladder5"},
        {"barcode": wl[21], "name": "Ladder6"},
        {"barcode": wl[22], "name": "Ladder7"},
        {"barcode": wl[23], "name": "Ladder8"},
    ]
).to_csv("bc.csv", sep="\t", index=False, header=None)
# %%
