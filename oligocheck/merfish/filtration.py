import numpy as np
import polars as pl

from oligocheck.algorithms import find_overlap, find_overlap_weighted


def handle_overlap(
    df: pl.DataFrame,
    criteria: list[pl.Expr],
    overlap: int = -2,
    n: int = 100,
):
    if len(gene := df.select(pl.col("gene").unique())) > 1:
        raise ValueError("More than one gene in filtered")
    gene = gene.item()

    df = df.sort(by=["is_ori_seq", "transcript_ori", "pos_end", "tm"], descending=[True, False, False, True])

    if not criteria:
        criteria = [pl.col("*")]

    ddf = df.lazy().with_row_count("index").with_columns(priority=pl.lit(0, dtype=pl.UInt8))
    priority = len(criteria)
    for criterion in reversed(criteria):
        ddf = ddf.update(
            ddf.filter(criterion).with_columns(priority=pl.lit(priority, dtype=pl.UInt8)), on="index"
        )
        priority -= 1
    ddf = ddf.collect()
    ddf.write_parquet(f"output/{gene}_filtered.parquet")
    df = ddf.filter(pl.col("priority") > 0)

    selected_global = set()
    tss = df["transcript_ori"].unique().to_list()
    for ts in tss:
        this_transcript = df.filter(pl.col("transcript_ori") == ts)
        print(len(this_transcript))
        if not len(this_transcript):
            continue

        for i in range(1, len(criteria) + 1):
            run = (
                df.filter((pl.col("priority") <= i) & ~pl.col("index").is_in(selected_global))
                .select(
                    [
                        "index",
                        "pos_start",
                        "pos_end",
                        "priority",
                    ]
                )
                .sort(["pos_end", "pos_start"])
            )
            if not len(run):
                continue

            priorities = np.sqrt(len(criteria) + 1 - run["priority"])
            try:
                if i == 1:
                    ols = find_overlap(run["pos_start"], run["pos_end"], overlap=overlap)
                else:
                    ols = find_overlap_weighted(run["pos_start"], run["pos_end"], priorities, overlap=overlap)
                sel_local = set(run[ols]["index"].to_list())
                print(i, len(sel_local))
                if len(sel_local) > n:
                    break
            except RecursionError:
                print("Recursion error")
                break

        selected_global |= sel_local  # type: ignore

    return df.filter(pl.col("index").is_in(selected_global))


def the_filter(df: pl.DataFrame, overlap: int = -1) -> pl.DataFrame:
    out = []
    for _, group in df.groupby("gene"):
        out.append(
            handle_overlap(
                group,
                criteria=[
                    # fmt: off
                    (pl.col("oks") > 4) & (pl.col("hp") < 35) & pl.col("match_max_all").lt(21) & pl.col("length").lt(42) & (pl.col("maps_to_pseudo").is_null() | pl.col("maps_to_pseudo").eq("")),
                    (pl.col("oks") > 4) & (pl.col("hp") < 35) & pl.col("match_max_all").lt(21) & pl.col("length").lt(42),
                    (pl.col("oks") > 3) & (pl.col("hp") < 35) & pl.col("match_max_all").lt(24),
                    (pl.col("oks") > 3) & pl.col("match_max_all").lt(28),
                    (pl.col("oks") > 2) & pl.col("match_max_all").lt(28),
                    (pl.col("oks") > 3),
                    # (pl.col("oks") > 1),
                    # fmt: on
                ],
                overlap=overlap,
            ).sort("priority")
        )
    return pl.concat(out)


def check_kmers(seq: str, kmer: set[str], n: int) -> bool:
    for i in range(len(seq) - n + 1):
        if seq[i : i + n] in kmer:
            return True
    return False
