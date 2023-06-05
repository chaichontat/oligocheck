# %%
import logging
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable

import click
import numpy as np
import polars as pl

from oligocheck.algorithms import find_overlap, find_overlap_weighted
from oligocheck.merfish.alignment import gen_fasta, gen_fastq, run_bowtie
from oligocheck.merfish.crawler import crawler
from oligocheck.merfish.external_data import ExternalData, get_rrna
from oligocheck.merfish.filtration import the_filter
from oligocheck.seqcalc import hp_fish, tm_fish
from oligocheck.sequtils import parse_sam, tm_hybrid

sys.setrecursionlimit(5000)
pl.Config.set_fmt_str_lengths(100)
pl.Config.set_tbl_rows(100)

fpkm = pl.read_parquet("data/fpkm/P0_combi.parquet")

try:
    profile
except NameError:
    profile = lambda x: x


# %%
def count_genes(df: pl.DataFrame) -> pl.DataFrame:
    return df.groupby("gene").count().sort("count")


def filter_gene(df: pl.DataFrame, gene: str) -> pl.DataFrame:
    return df.filter(pl.col("gene") == gene)


# GENCODE primary only
gtf = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts.parquet",
    path="data/mm39/gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf",
    fasta="data/mm39/combi.fa.gz",
)

gtf_all = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts_all.parquet",
    path="data/mm39/Mus_musculus.GRCm39.109.gtf",
    fasta="data/mm39/combi.fa.gz",
)

rrnas = get_rrna("data/mm39/Mus_musculus.GRCm39.ncrna.fa.gz")
trna_rna_kmers = set(
    pl.read_csv("data/mm39/trcombi.txt", separator=" ", has_header=False, new_columns=["kmer", "count"])[
        "kmer"
    ]
)


@profile
def block_bowtie(gene: str, tss: Iterable[str], temp: Path):
    out = []
    for ts in tss:
        seq = gtf.get_seq(ts)
        res = crawler(seq, prefix=f"{gene}_{ts}")
        out.append(
            run_bowtie(
                gen_fastq(res["name"], res["seq"]).getvalue(),
                "data/mm39/mm39",
                seed_length=15,
                threshold=17,
                n_return=-1,
            )
        )
    return out, res


def get_pseudogenes(gene: str, y: pl.DataFrame):
    counts = (
        y.groupby("transcript")
        .count()
        .join(fpkm, left_on="transcript", right_on="transcript_id(s)", how="left")
        .join(
            gtf_all[["transcript_id", "transcript_name"]],
            left_on="transcript",
            right_on="transcript_id",
        )
        .sort("count", descending=True)
    )

    # Filtered based on expression and number of probes aligned.
    ok = counts.filter(
        (pl.col("count") > 0.1 * pl.col("count").first())
        & (
            pl.col("FPKM") < 0.05 * pl.col("FPKM").first()
            if counts[0, "FPKM"] is not None and counts[0, "FPKM"] > 1.0
            else pl.lit(True)
        )
        & (pl.col("transcript_name").str.starts_with(gene) | pl.col("transcript_name").str.starts_with("Gm"))
    )
    return ok["transcript"]


def filter_match(y: pl.DataFrame, acceptable_tss: Iterable[str], match: float = 0.8, match_max: float = 0.6):
    nogo = y.filter(
        (pl.col("match").gt(pl.col("length") * match))
        & (pl.col("match_max").gt(pl.col("length") * match_max))
        & ~pl.col("transcript").str.contains("|".join([f"({x})" for x in acceptable_tss]))
    )
    return y.filter(~pl.col("name").is_in(nogo["name"]))


# %%


# %%
def run(gene: str, overlap: int = -2, allow_pseudo: bool = True):
    Path("temp").mkdir(exist_ok=True)
    if Path(f"output/{gene}_final.parquet").exists():
        return
    # wants = list(filter(lambda x: x, Path("celltypegenes.csv").read_text().splitlines()))
    # gene = wants[7]
    tss_gencode = gtf.filter_gene(gene)["transcript_id"]
    tss_all = gtf_all.filter_gene(gene)["transcript_id"]
    longest = max(tss_gencode, key=lambda x: len(gtf.get_seq(x)))

    if not Path(f"output/{gene}_all.parquet").exists():
        bt, crawled = block_bowtie(gene, [longest], Path("temp"))
        y = parse_sam(bt[0])
        offtargets = (
            y["transcript"]
            .value_counts()
            .sort("counts", descending=True)
            .with_columns(name=pl.col("transcript").apply(gtf_all.ts_to_gene))
        )

        y: pl.DataFrame = (
            y.join(
                y[["id", "mismatched_reference"]]
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
            ).with_columns(revcomped=pl.col("flag").map(lambda col: np.bitwise_and(col, 16) > 0))
            # .filter(~pl.col("revcomped"))
        )

        y.write_parquet(f"output/{gene}_all.parquet")
        offtargets.write_parquet(f"output/{gene}_offtargets.parquet")
    else:
        print(f"{gene} reading from cache.")
        y = pl.read_parquet(f"output/{gene}_all.parquet")
        offtargets = pl.read_parquet(f"output/{gene}_offtargets.parquet")

    tss_pseudo = get_pseudogenes(gene, y) if allow_pseudo else pl.Series()
    ff = filter_match(y, [*tss_all, *tss_pseudo], match=0.8, match_max=0.6)

    if len(tss_pseudo):
        names_with_pseudo = ff.filter(pl.col("transcript").str.contains("|".join(tss_pseudo))).rename(
            dict(transcript="maps_to_pseudo")
        )[["name", "maps_to_pseudo"]]
    else:
        names_with_pseudo = pl.DataFrame({"name": ff["name"].unique(), "maps_to_pseudo": pl.lit(False)})

    ff = (
        ff.lazy()
        .filter("is_ori_seq")
        .with_columns(
            [
                pl.col("transcript").apply(gtf_all.ts_to_gene).alias("transcript_name"),
                pl.col("seq").str.contains("GGGG").is_not().alias("ok_quad_c"),
                pl.col("seq").str.contains("TTTT").is_not().alias("ok_quad_a"),
                (pl.col("seq").str.count_match("T") / pl.col("seq").str.n_chars() < 0.28).alias("ok_comp_a"),
                pl.all(
                    [pl.col("seq").str.slice(-6 - i, 6).str.count_match("G").lt(4) for i in range(6)]
                ).alias("ok_stack_c"),
                (pl.col("seq").str.count_match("G|C") / (pl.col("seq").str.n_chars())).alias("gc_content"),
            ]
        )
        .with_columns(
            [
                (pl.col("gc_content").is_between(0.35, 0.65)).alias("ok_gc"),
                pl.col("seq").apply(tm_hybrid).alias("tm"),
                pl.col("seq").apply(hp_fish).alias("hp") - 0.65 * 30,
            ]
        )
        .with_columns(oks=pl.sum(pl.col("^ok_.*$")))
        .collect()
        .join(names_with_pseudo, on="name", how="left")
    )
    ff.unique("name").write_parquet(f"output/{gene}_filtered.parquet")
    final = the_filter(ff, gtf, overlap=overlap)
    assert not (final["flag"].to_numpy() & 16).any()
    final.write_parquet(f"output/{gene}_final.parquet")


@click.command()
@click.argument("gene")
# @click.option("--output", "-o", type=click.Path(), default="output/")
@click.option("--debug", "-d", is_flag=True)
# @click.option("--allow-pseudo", "-p", is_flag=True)
@click.option("--overlap", "-O", type=int, default=-2)
def main(gene: str, overlap: int = -2, allow_pseudo: bool = True, debug: bool = False):
    if debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        run(gene, overlap=overlap, allow_pseudo=allow_pseudo)
    except Exception as e:
        raise Exception(f"Failed to run {gene}") from e

    # Path(output).mkdir(exist_ok=True, parents=True)
    # m.write_parquet(f"{output}/{gene}.parquet")
    # with open(f"{output}/{gene}.json", "w") as f:
    #     json.dump(md, f)


if __name__ == "__main__":
    main()
