# %%
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from oligocheck.merfish.external_data import gene_to_eid, get_gencode

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
    tss = tuple(everything[everything.gene_id == eid].index)
    if not criteria:
        criteria = [lambda _: True]

    selected = {}
    for ts in tss:
        this_transcript = filtered[filtered["transcript_ori"] == ts]
        if len(this_transcript) == 0:
            print("No match found for", ts)
            continue
        forbidden = np.zeros(this_transcript.pos_end.max() + 1 + max(0, overlap), dtype=bool)
        priority = 1
        for criterion in criteria:
            if len(selected) >= n:
                break
            for idx, r in this_transcript.iterrows():
                if idx in selected or not criterion(r):
                    continue
                if np.any(forbidden[r.pos_start - 1 : r.pos_end + 1 - overlap]):
                    continue
                selected[idx] = priority
                forbidden[r.pos_start - 1 : r.pos_end + 1 - overlap] = 1
            priority += 1
    filtered = filtered.loc[selected.keys()]
    filtered["priority"] = filtered.index.map(selected)
    return filtered


def the_filter(df: pd.DataFrame, genes: Iterable[str], overlap: int = -1):
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
                    & (x.nonspecific_binding < 0.05),
                    lambda x: (x.tm > 49)
                    & (x.tm < 54)
                    & (x.oks > 3)
                    & (x.hp < 35)
                    & (x.nonspecific_binding < 0.001),
                    lambda x: (x.tm > 47) & (x.tm < 56) & (x.hp < 40) & (x.oks > 3),
                    lambda x: (x.tm > 46) & (x.tm < 56) & (x.hp < 40) & (x.oks > 3),
                    lambda x: (x.tm > 46) & (x.tm < 56) & (x.hp < 40) & (x.oks > 2),
                ],
                overlap=overlap,
                n=100,
            )
        )
    return pd.concat(out)
