from typing import Iterable

import polars as pl


class GeneFrame(pl.DataFrame):
    def __init__(self, df: pl.DataFrame):
        self._df = df._df

    def gene(self, gene: str):
        return self.filter(pl.col("gene") == gene)

    def count(self, col: str = "gene"):
        return self.groupby(col).agg(pl.count()).sort("count")

    def filter_(self, col: str, value: str):
        return self.filter(pl.col(col) == value)

    @classmethod
    def concat(cls, dfs: Iterable[pl.DataFrame]):
        return cls(pl.concat(dfs))

    @classmethod
    def read_parquet(cls, path: str):
        return cls(pl.read_parquet(path))
