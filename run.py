#%%
from itertools import combinations
from pathlib import Path

from oligocheck.io import get_whitelist
from oligocheck.tenx import gen_oligo
from oligocheck.utils import hamming

whitelist = get_whitelist()
bcs = list(range(1, 11))
for x, y in combinations(bcs, 2):
    assert hamming(whitelist[x - 1], whitelist[y - 1]) > 4, f"{x} {y} distance less than 5."

Path("out.txt").write_text(gen_oligo(bcs, "polyA", "ADT"))

# %%
