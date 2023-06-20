# %%

import polars as pl

from oligocheck.seqcalc import hp_fish, tm_hybrid
from oligocheck.sequtils import gc_content


def crawler(
    seq: str,
    prefix: str,
    length_limit: tuple[int, int] = (25, 46),
    gc_limit: tuple[float, float] = (0.3, 0.7),
    tm_limit: float = 51,
    hairpin_limit: float = 40 + 0.65 * 30,
) -> pl.DataFrame:
    j = length_limit[0]
    names = []
    seqs = []

    for i in range(len(seq) - length_limit[0]):
        while True:
            if j - i < length_limit[0]:
                j += 1
                continue
            if j - i > length_limit[1] or j > len(seq):
                break
            if (gc_content(seq[i:j]) < gc_limit[0]) or (gc_content(seq[i:j]) > gc_limit[1]):
                j += 1
                continue

            if tm_hybrid(seq[i:j]) > tm_limit:
                if hp_fish(seq[i:j]) < hairpin_limit:
                    names.append(f"{prefix}:{i}-{j-1}")
                    seqs.append(seq[i:j])
                break
            j += 1
    return pl.DataFrame(dict(name=names, seq=seqs)).filter(
        ~pl.col("seq").str.contains(r"AAAAA|TTTTT|CCCCC|GGGGG")
    )


# %%
# %%
# %%
# %%
# %%
