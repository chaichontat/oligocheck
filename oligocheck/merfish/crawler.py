# %%

import polars as pl
import primer3

from oligocheck.sequtils import gc_content, tm_hybrid


def crawler(seq: str, prefix: str) -> pl.DataFrame:
    minmax = [24, 41]
    gc_limit = [0.3, 0.7]
    j = minmax[0]
    hp_limit = 40 + 0.65 * 30

    names = []
    seqs = []

    for i in range(len(seq) - minmax[0]):
        while True:
            if j - i < minmax[0]:
                j += 1
                continue
            if j - i > minmax[1] or j > len(seq):
                break
            if (gc_content(seq[i:j]) < gc_limit[0]) or (gc_content(seq[i:j]) > gc_limit[1]):
                j += 1
                continue

            if (tm := tm_hybrid(seq[i:j])) > 49:
                if (
                    primer3.calc_hairpin_tm(seq[i:j], mv_conc=330, dv_conc=0, dntp_conc=0, dna_conc=1)
                    < hp_limit
                ):
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
