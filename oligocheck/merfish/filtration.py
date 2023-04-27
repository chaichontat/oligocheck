import numpy as np
import polars as pl

from oligocheck.merfish.external_data import ExternalData


def handle_overlap(
    ensembl: ExternalData,
    df: pl.DataFrame,
    criteria: list[pl.Expr],
    overlap: int = -1,
    n: int = 200,
):
    if len(df.select(pl.col("gene").unique())) > 1:
        raise ValueError("More than one gene in filtered")
    gene = df["gene"].unique().item()
    df = df.sort(
        by=["is_ori_seq", "transcript_ori", "pos_start", "tm"], descending=[True, False, False, True]
    )
    eid = ensembl.gene_to_eid(gene)
    tss = tuple(ensembl.filter(pl.col("gene_id") == eid)["transcript_id"])

    if not criteria:
        criteria = [pl.all("*")]

    df = df.with_columns(
        [
            pl.lit(0, dtype=pl.UInt8).alias("priority"),
            pl.arange(0, len(df), dtype=pl.UInt32).alias("index"),
        ]
    )
    selected = set()

    for ts in tss:
        this_transcript = df.filter(pl.col("transcript_ori") == ts)
        if len(this_transcript) == 0:
            print("No match found for", ts)
            continue

        forbidden = np.zeros(df.select(pl.col("pos_end")).max().item() + 1 + max(0, overlap), dtype=bool)
        priority = 1
        for criterion in criteria:
            if len(selected) >= n:
                break
            for idx, st, end in (
                df.filter(criterion & ~pl.col("index").is_in(selected))
                .select(["index", "pos_start", "pos_end"])
                .iter_rows()
            ):
                if np.any(forbidden[st - 1 : end + 1 - overlap]):
                    continue
                selected.add(idx)
                df[idx, "priority"] = priority
                forbidden[st - 1 : end + 1 - overlap] = 1
            priority += 1

    return df.filter(pl.col("priority") > 0)
