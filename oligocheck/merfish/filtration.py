import polars as pl

from oligocheck.algorithms import find_overlap, find_overlap_weighted
from oligocheck.merfish.external_data import ExternalData


def count_match(df: pl.DataFrame) -> pl.DataFrame:
    return df.join(
        df[["id", "mismatched_reference"]]
        .with_columns(mismatched_reference=pl.col("mismatched_reference").str.extract_all(r"(\d+)"))
        .explode("mismatched_reference")
        .with_columns(pl.col("mismatched_reference").cast(pl.UInt8))
        .groupby("id")
        .agg(
            match=pl.col("mismatched_reference").sum(),
            match_max=pl.col("mismatched_reference").max(),
        ),
        on="id",
        how="left",
    )


def handle_overlap(
    ensembl: ExternalData,
    df: pl.DataFrame,
    criteria: list[pl.Expr],
    overlap: int = -1,
    n: int = 100,
):
    if len(gene := df.select(pl.col("gene").unique())) > 1:
        raise ValueError("More than one gene in filtered")
    gene = gene.item()
    df = df.sort(by=["is_ori_seq", "transcript_ori", "pos_end", "tm"], descending=[True, False, False, True])
    eid = ensembl.gene_to_eid(gene)
    tss = tuple(ensembl.filter(pl.col("gene_id") == eid)["transcript_id"])

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

            priorities = (5 - run["priority"]).to_list()
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


def the_filter(df: pl.DataFrame, gtf: ExternalData, overlap: int = -1) -> pl.DataFrame:
    out = []
    for _, group in df.groupby("gene"):
        out.append(
            handle_overlap(
                gtf,
                group,
                criteria=[
                    # fmt: off
                    (pl.col("oks") > 4) & (pl.col("hp") < 35)& pl.col("tm").is_between(50, 54),
                    (pl.col("oks") > 3) & (pl.col("hp") < 35) & pl.col("tm").is_between(50, 54),
                    (pl.col("oks") > 3) & (pl.col("hp") < 40) & pl.col("tm").is_between(49, 56),
                    (pl.col("oks") > 2) & (pl.col("hp") < 40) & pl.col("tm").is_between(49, 56),
                    # fmt: on
                ],
                overlap=overlap,
            )
        )
    return pl.concat(out)
