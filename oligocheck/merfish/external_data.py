import fcntl
from contextlib import contextmanager
from functools import cache
from pathlib import Path
from typing import overload

import mygene
import pandas as pd
import polars as pl
import pyfastx
from Bio import SeqIO

mg = mygene.MyGeneInfo()


@contextmanager
def lock(path: str):
    locked_file_descriptor = open(Path(path), "w+")
    fcntl.lockf(locked_file_descriptor, fcntl.LOCK_EX)
    yield locked_file_descriptor
    locked_file_descriptor.close()


class ExternalData:
    def __init__(self, *, cache: Path | str, path: Path | str, fasta: Path | str) -> None:
        self.fa = pyfastx.Fasta(fasta, key_func=lambda x: x.split(" ")[0].split(".")[0])

        if Path(cache).exists():
            self.gtf: pl.DataFrame = pl.read_parquet(cache)
            return
        self.gtf = self.parse_gtf(path)
        self.gtf.write_parquet(cache)

    @cache
    def gene_info(self, gene: str) -> pl.DataFrame:
        return self.gtf.filter(pl.col("gene_name") == gene)

    @cache
    def gene_to_eid(self, gene: str) -> str:
        return self.gene_info(gene)[0, "gene_id"]

    @cache
    def eid_to_ts(self, eid: str) -> str:
        eid = eid.split(".")[0]
        return self.gtf.filter(pl.col("transcript_id") == eid)[0, "transcript_id"]

    @cache
    def all_transcripts(self, gene: str | None = None, *, eid: str | None = None) -> pl.Series:
        if gene is not None:
            return self.gtf.filter(pl.col("gene_name") == gene)["transcript_id"]
        return self.gtf.filter(pl.col("gene_id") == eid)["transcript_id"]

    @cache
    def get_seq(self, eid: str) -> str:
        res = self.fa[eid.split(".")[0]].seq
        if not res:
            raise ValueError(f"Could not find {eid}")
        return res

    @overload
    def __getitem__(self, eid: str) -> pl.Series:
        ...

    @overload
    def __getitem__(self, eid: list[str]) -> pl.DataFrame:
        ...

    def __getitem__(self, eid: str | list[str]):
        return self.gtf[eid]

    @staticmethod
    def parse_gtf(path: str | Path) -> pl.DataFrame:
        # fmt: off
        # To get the original keys.
        # list(reduce(lambda x, y: x | json.loads(y), jsoned['jsoned'].to_list(), {}).keys())
        attr_keys = (
            "gene_id", "transcript_id", "gene_type", "gene_name", "transcript_type",
            "transcript_name", "level", "transcript_support_level", "mgi_id", "tag",
            # "havana_gene", "havana_transcript", "protein_id", "ccdsid", "ont",
        )
        # fmt: on

        return (
            pl.scan_csv(
                path,
                comment_char="#",
                separator="\t",
                has_header=False,
                new_columns=[
                    "seqname",
                    "source",
                    "feature",
                    "start",
                    "end",
                    "score",
                    "strand",
                    "frame",
                    "attribute",
                ],
                dtypes=[pl.Utf8, pl.Utf8, pl.Utf8, pl.UInt32, pl.UInt32, pl.Utf8, pl.Utf8, pl.Utf8, pl.Utf8],
            )
            .filter(pl.col("feature") == "transcript")
            .with_columns(
                pl.concat_str(
                    [
                        pl.lit("{"),
                        pl.col("attribute")
                        .str.replace_all(r"; (\w+) ", r', "$1": ')
                        .str.replace_all(";", ",")
                        .str.replace(r"(\w+) ", r'"$1": ')
                        .str.replace(r",$", ""),
                        pl.lit("}"),
                    ]
                ).alias("jsoned")
            )
            .with_columns(
                [
                    pl.col("jsoned").str.json_path_match(f"$.{name}").alias(name)
                    # .cast(pl.Categorical if "type" in name or "tag" == name else pl.Utf8)
                    for name in attr_keys
                ]
            )
            # .drop(["attribute", "jsoned"])
            .with_columns(
                [
                    pl.col("gene_id").str.extract(r"(\w+)(\.\d+)?").alias("gene_id"),
                    pl.col("transcript_id").str.extract(r"(\w+)(\.\d+)?").alias("transcript_id"),
                ]
            )
            .collect()
        )


def get_rrna(path: str) -> set[str]:
    if not Path(path).exists():
        # Get tRNA and rRNA
        out = []
        ncrna = SeqIO.parse("data/mm39/Mus_musculus.GRCm39.ncrna.fa", "fasta")
        for line in ncrna:
            attrs = line.description.split(", ")[0].split(" ")
            actual = attrs[:7]
            description = " ".join(attrs[7:])
            attrs = [":".join(x.split(":")[1:]) if ":" in x else x for x in actual] + [
                str(line.seq),
                description,
            ]
            if len(attrs) < 8:
                attrs.append("")
            out.append(attrs)

        out = pd.DataFrame.from_records(
            out,
            columns=[
                "transcript_id",
                "type",
                "pos",
                "gene_id",
                "gene_biotype",
                "transcript_biotype",
                "gene_symbol",
                "seq",
                "description",
            ],
        )
        with open(path, "w") as f:
            for _, row in out[out.gene_biotype == "rRNA"].iterrows():
                f.write(f">{row.transcript_id}\n{row.seq}\n")

    return set(pd.read_table(path, header=None)[0].str.split(">", expand=True)[1].dropna())
