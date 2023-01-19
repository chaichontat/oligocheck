#%%
from itertools import combinations
from pathlib import Path

from oligocheck.io import get_whitelist
from oligocheck.tenx import gen_oligo
from oligocheck.utils import hamming

whitelist = get_whitelist()
bcs = list(range(1001, 1025))
for x, y in combinations(bcs, 2):
    assert hamming(whitelist[x - 1], whitelist[y - 1]) > 4, f"{x} {y} distance less than 5."


def builder(seq: str, i: int) -> str:
    return seq.replace("*", "")


def namer(seq: str, i: int) -> str:
    return f"Hash-sci-{i}-nextera"


Path("out.txt").write_text(gen_oligo(bcs, builder=builder, namer=namer, provider="sci"))

# %%
