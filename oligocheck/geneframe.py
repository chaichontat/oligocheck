from functools import reduce
from io import StringIO
from typing import Collection, Iterable, TypeVar

import polars as pl

from oligocheck.seqcalc import tm_match

T = TypeVar("T", bound=pl.DataFrame)


class GeneFrame(pl.DataFrame):
    # NECESSARY_COLUMNS = {"seq", "transcript", "pos_start", "pos_end"}

    def __init__(self, df: pl.DataFrame):
        self._df = df._df

    def gene(self, gene: str):
        return self.filter(pl.col("gene") == gene)

    def count(self, col: str = "gene", descending: bool = False):
        return self.groupby(col).agg(pl.count()).sort("count", descending=descending)

    def filter_eq(self, **kwargs: str | float):
        return self.filter(reduce(lambda x, y: x & y, [pl.col(k) == v for k, v in kwargs.items()]))

    def filter_isin(self, **kwargs: Collection[str] | pl.Series):
        return self.filter(reduce(lambda x, y: x & y, [pl.col(k).is_in(v) for k, v in kwargs.items()]))

    def left_join(
        self,
        other: pl.DataFrame,
        on: str | None = None,
        left_on: str | None = None,
        right_on: str | None = None,
    ):
        return self.join(other, on=on, left_on=left_on, right_on=right_on, how="left")

    @classmethod
    def concat(cls, dfs: Iterable[pl.DataFrame]):
        return cls(pl.concat(dfs))

    @classmethod
    def read_parquet(cls, path: str):
        return cls(pl.read_parquet(path))

    @staticmethod
    def _count_match(df: T) -> T:
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

    @classmethod
    def from_sam(cls, sam: str, split_name: bool = True, count_match: bool = True):
        # s = (
        #     pl.DataFrame(dict(strs=[sam]))
        #     .lazy()
        #     .with_columns(pl.col("strs").str.split("\n"))
        #     .explode("strs")
        #     .with_columns(pl.col("strs").str.strip().str.split("\t").arr.slice(0, 10))
        #     .with_row_count("id")
        #     .explode("strs")
        #     .with_columns(col_nm="string_" + pl.arange(0, pl.count()).cast(pl.Utf8).str.zfill(2).over("id"))
        #     .sort(["col_nm", "id"])
        #     .collect()
        # )

        # faster than pivot
        # key_optional = {
        #     "AS": "aln_score",
        #     "XS": "aln_score_best",
        #     "XN": "n_ambiguous",
        #     "XM": "n_mismatches",
        #     "XO": "n_opens",
        #     "XG": "n_extensions",
        #     "NM": "edit_distance",
        #     "MD": "mismatched_reference",
        #     "YT": "pair_state",
        # }

        df = (
            pl.read_csv(StringIO(sam), separator="\n", has_header=False)
            .lazy()
            .with_row_count("id")
            .with_columns(temp=pl.col("column_1").str.split_exact("\t", 9))
            .unnest("temp")
            .rename(
                {
                    f"field_{i}": x
                    for i, x in enumerate(
                        [
                            "name",
                            "flag",
                            "transcript",
                            "pos",
                            "mapq",
                            "cigar",
                            "rnext",
                            "pnext",
                            "tlen",
                            "seq",
                        ]
                    )
                }
            )
            .with_columns(
                flag=pl.col("flag").cast(pl.UInt16),
                pos=pl.col("pos").cast(pl.UInt32),
                mapq=pl.col("mapq").cast(pl.UInt8),
                aln_score=pl.col("column_1").str.extract(r"AS:i:(\d+)").cast(pl.UInt16),
                aln_score_best=pl.col("column_1").str.extract(r"XS:i:(\d+)").cast(pl.UInt16),
                n_ambiguous=pl.col("column_1").str.extract(r"XN:i:(\d+)").cast(pl.UInt16),
                n_mismatches=pl.col("column_1").str.extract(r"XM:i:(\d+)").cast(pl.UInt16),
                n_opens=pl.col("column_1").str.extract(r"XO:i:(\d+)").cast(pl.UInt16),
                n_extensions=pl.col("column_1").str.extract(r"XG:i:(\d+)").cast(pl.UInt16),
                edit_distance=pl.col("column_1").str.extract(r"NM:i:(\d+)").cast(pl.UInt16),
                mismatched_reference=pl.col("column_1").str.extract(r"MD:Z:(\S+)"),
            )
            .drop(["column_1", "mapq", "rnext", "pnext", "tlen"])
            .with_columns(
                [
                    pl.when(pl.col("transcript").str.contains(r"(.*)\.\d+"))
                    .then(pl.col("transcript").str.extract(r"(.*)\.\d+"))
                    .otherwise(pl.col("transcript"))
                    .alias("transcript")
                ]
                + [
                    pl.col("name").str.extract(r"(.+)_(.+):(\d+)-(\d+)", 1).alias("gene"),
                    pl.col("name").str.extract(r"(.+)_(.+):(\d+)-(\d+)", 2).alias("transcript_ori"),
                    pl.col("name")
                    .str.extract(r"(.+)_(.+):(\d+)-(\d+)", 3)
                    .cast(pl.UInt32)
                    .alias("pos_start"),
                    pl.col("name").str.extract(r"(.+)_(.+):(\d+)-(\d+)", 4).cast(pl.UInt32).alias("pos_end"),
                ]
                if split_name
                else []
            )
            .with_columns(
                [
                    (pl.col("transcript") == pl.col("transcript_ori")).alias("is_ori_seq"),
                ]
                + [
                    (pl.col("pos_end") - pl.col("pos_start") + 1).alias("length"),
                ]
                if split_name
                else []
            )
            .collect()
        )
        return cls._count_match(cls(df)) if count_match else cls(df)

    def filter_by_match(
        self,
        acceptable_tss: Collection[str],
        match: float = 0.8,
        match_max: float = 0.7,
        return_nogo: bool = False,
    ):
        nogo = self.filter(
            (pl.col("match").gt(pl.col("length") * match))
            & (pl.col("match_max").gt(pl.col("length") * match_max))
            & ~pl.col("transcript").is_in(acceptable_tss)
        )

        tm_offtarget = self.filter(
            pl.col("match_max").is_between(pl.col("length") * 0.5, pl.col("length") * match_max + 0.01)
            & pl.col("match_max").gt(15)
        ).with_columns(
            tm_offtarget=pl.struct(["seq", "cigar", "mismatched_reference"]).apply(
                lambda x: tm_match(x["seq"], x["cigar"], x["mismatched_reference"])  # type: ignore
            )
        )
        unique = tm_offtarget.groupby("name").agg(
            pl.max("tm_offtarget").alias("max_tm_offtarget").cast(pl.Float32),
            pl.max("match_max").alias("match_max_all"),
        )
        nogo_soft = tm_offtarget.filter(pl.col("tm_offtarget").gt(40))

        # print(f"Found {len(nogo.unique('name'))} hard and {len(nogo_soft.unique('name'))} soft no-gos.")
        filtered = (
            self.filter(~pl.col("name").is_in(nogo["name"]) & ~pl.col("name").is_in(nogo_soft["name"]))
            .join(unique, on="name", how="left")
            .with_columns(
                max_tm_offtarget=pl.col("max_tm_offtarget").fill_null(0.0),
                match_max_all=pl.col("match_max_all").fill_null(0),
            )
        )
        return (
            filtered
            if not return_nogo
            else (filtered, nogo["name"].unique().to_list() + nogo_soft["name"].unique().to_list())
        )


def concat(dfs: Iterable[pl.DataFrame]) -> GeneFrame:
    return GeneFrame(pl.concat(dfs))
