# %%
from pathlib import Path

import numpy as np
import polars as pl
from expression import pipe
from Levenshtein import distance

from oligocheck.merfish.alignment import gen_fastq, run_bowtie
from oligocheck.sequtils import parse_sam

seed = "TTACACTCCATCCACTCAA"

for f in Path("generator/").glob("combinations_*.txt"):
    print(f)
    f.with_suffix(".sam").write_bytes(
        run_bowtie(
            f.read_bytes(),
            "data/humouse/humouse",
            seed_length=13,
            n_return=1,
            threshold=15,
        )
    )
# %%

sam = pl.concat(parse_sam(f.read_text()) for f in Path("generator/").glob("combinations_*.sam"))
# %%
ok = sam.filter(pl.col("transcript") == "*")
ok.write_parquet("generator/lessthan30.parquet")
ori = ok.select(name="name", column_1="seq")

# %%


def run_filter(
    df: pl.DataFrame,
    *,
    seed_length: int,
    n_return: int,
    threshold: int,
) -> pl.DataFrame:
    if "aln_score" in df.columns:
        # grouped = df.groupby("name").agg(pl.max("aln_score"))
        d = df.join(ori, on="name", how="inner")
    else:
        d = df
    print(len(d))
    res = pipe(
        gen_fastq(d["name"], d["column_1"]),
        lambda x: run_bowtie(
            x.getvalue(),
            "data/humouse/humouse",
            seed_length=seed_length,
            n_return=n_return,
            threshold=threshold,
        ),
        lambda x: parse_sam(x, split_name=False),
    )
    return res


# %%
def histogram(df: pl.DataFrame) -> None:
    df.groupby("name").agg(pl.max("aln_score")).to_pandas().hist()


def get_unaligned(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter(pl.col("transcript") == "*")


# def filter_seed(df: pl.DataFrame, threshold: int = 8) -> pl.DataFrame:
#     d = ori[df["name"].unique().cast(pl.UInt32)].with_columns(
#         dist=pl.col("column_1").apply(lambda x: distance(seed, x))
#     )
#     return d.filter(pl.col("dist") > threshold)


def filter_max_score(df: pl.DataFrame, threshold: int):
    return df.groupby("name").agg(pl.max("aln_score")).filter(pl.col("aln_score") <= threshold)


# %%
res = run_filter(ori, seed_length=12, n_return=5, threshold=16)
# %%
res2 = run_filter(get_unaligned(res), seed_length=11, n_return=10, threshold=16)
# %%
res3 = run_filter(get_unaligned(res2), seed_length=10, n_return=50, threshold=16)

# %%
res4 = run_filter(get_unaligned(res3), seed_length=10, n_return=100, threshold=14)


# %%
res5 = run_filter(
    res4.groupby("name").agg(pl.max("aln_score")).filter(pl.col("aln_score") <= 31),
    seed_length=10,
    n_return=500,
    threshold=15,
)
# %%

res6 = run_filter(
    res5.groupby("name")
    .agg(pl.max("aln_score"))
    .filter(pl.col("aln_score").lt(32) | pl.col("aln_score").is_null()),
    seed_length=9,
    n_return=100,
    threshold=15,
)

# %%%
filtered = ori.join(
    (
        res6.groupby("name")
        .agg(pl.max("aln_score"))
        .filter(pl.col("aln_score").lt(32) | pl.col("aln_score").is_null())
    ),
    on="name",
    how="inner",
)
filtered.write_parquet("generator/filtered.parquet")
# %%


def screen_seed(seed: int, dist: int = 7):
    avail = filtered["column_1"]
    selected = [avail[seed]]
    for i in range(len(avail)):
        if i == seed:
            continue

        if distance(avail[i], selected[-1]) < dist:
            continue

        # Distance check
        for s in selected:
            if distance(s, avail[i]) < dist:
                break
        else:
            selected.append(avail[i])
    return selected


y = [screen_seed(i) for i in range(25)]
best_seed = int(np.array(list(map(len, y))).argmax())
# %%
prefinal = screen_seed(best_seed)
prefinal = parse_sam(
    run_bowtie(
        gen_fastq(pl.Series(prefinal), pl.Series(prefinal)).getvalue(),
        "data/humouse/humouse",
        seed_length=9,
        n_return=-1,
        threshold=14,
    ),
    split_name=False,
)
histogram(prefinal)
# %%
final = (
    prefinal.groupby("name")
    .agg(pl.max("aln_score"))
    .filter(pl.col("aln_score") <= 31)
    .select(seq="name")
    .with_row_count("id", offset=1)
)
final.write_csv("generator/readouts.csv")
# %%
