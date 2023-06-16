# %%
import logging
import sys
from pathlib import Path
from typing import cast

import click
import numpy as np
import polars as pl

from oligocheck.geneframe import GeneFrame
from oligocheck.logging import setup_logging
from oligocheck.merfish.alignment import gen_fastq, run_bowtie
from oligocheck.merfish.crawler import crawler
from oligocheck.external.external_data import ExternalData
from oligocheck.merfish.filtration import check_kmers, the_filter
from oligocheck.seqcalc import hp_fish, tm_hybrid

setup_logging()

sys.setrecursionlimit(5000)
pl.Config.set_fmt_str_lengths(100)
pl.Config.set_tbl_rows(100)

fpkm = pl.read_parquet("data/fpkm/P0_combi.parquet")
kmer18 = pl.read_csv(
    "data/mm39/kmer_genome18.txt", separator=" ", has_header=False, new_columns=["kmer", "count"]
)
trna_rna_kmers = set(
    pl.read_csv(
        "data/mm39/kmer_trcombi15.txt", separator=" ", has_header=False, new_columns=["kmer", "count"]
    )["kmer"]
)
kmerset = set(kmer18["kmer"])

try:
    profile
except NameError:
    profile = lambda x: x


# %%
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


@profile
def block_bowtie(gene: str, ts: str):
    seq = gtf.get_seq(ts)
    crawled = crawler(seq, prefix=f"{gene}_{ts}")
    return (
        GeneFrame.from_sam(
            run_bowtie(
                gen_fastq(crawled["name"], crawled["seq"]).getvalue(),
                "data/mm39/mm39",
                seed_length=13,
                threshold=18,
                n_return=500,
            )
        ),
        crawled,
    )


def get_pseudogenes(gene: str, y: GeneFrame, limit: int = 5) -> tuple[pl.Series, pl.Series]:
    print(type(y.count("transcript")))
    counts = (
        y.count("transcript")
        .left_join(fpkm, left_on="transcript", right_on="transcript_id(s)")
        .join(
            gtf_all[["transcript_id", "transcript_name"]],
            left_on="transcript",
            right_on="transcript_id",
        )
        .sort("count", descending=True)
    )

    # Filtered based on expression and number of probes aligned.
    ok = counts.filter(
        (pl.col("count") > 0.1 * pl.col("count").max())
        & (pl.col("FPKM").lt(0.05 * pl.col("FPKM").first()) | pl.col("FPKM").lt(1) | pl.col("FPKM").is_null())
        & (
            pl.col("transcript_name").str.starts_with(gene + "-ps")
            | pl.col("transcript_name").str.starts_with("Gm")
        )
    )
    return ok[:limit]["transcript"], ok[:limit]["transcript_name"]


# %%


# %%
def run(
    gene: str,
    output: str | Path,
    overlap: int = -2,
    allow_pseudo: bool = True,
    ignore_revcomp: bool = False,
    realign: bool = False,
):
    (output := Path(output)).mkdir(exist_ok=True)
    tss_gencode = set(gtf.filter_gene(gene)["transcript_id"])
    tss_all = gtf_all.filter_gene(gene)["transcript_id"]
    if not len(tss_gencode):
        raise ValueError(f"Gene {gene} not found in GENCODE.")
    canonical = gtf_all.filter_gene(gene).filter(pl.col("tag") == "Ensembl_canonical")[0, "transcript_id"]
    tss_gencode.add(canonical)  # Some canonical transcripts are not in GENCODE.

    if realign or not (output / f"{gene}_all.parquet").exists():
        y, crawled = block_bowtie(gene, canonical)

        offtargets = (
            y["transcript"]
            .value_counts()
            .sort("counts", descending=True)
            .with_columns(name=pl.col("transcript").apply(gtf_all.ts_to_gene))
        )

        y.write_parquet(output / f"{gene}_all.parquet")
        offtargets.write_parquet(output / f"{gene}_offtargets.parquet")
    else:
        print(f"{gene} reading from cache.")
        y = GeneFrame.read_parquet(output / f"{gene}_all.parquet")
        offtargets = pl.read_parquet(output / f"{gene}_offtargets.parquet")

    if ignore_revcomp:
        print("Ignoring revcomped reads.", len(y), end=" ")
        y = y.filter(pl.col("flag").map(lambda col: np.bitwise_and(col, 16) == 0))
        print(len(y))

    # Print most common offtargets
    print(
        y.count("transcript", descending=True)
        .filter(pl.col("count") > 0.1 * pl.col("count").first())
        .with_columns(name=pl.col("transcript").apply(gtf_all.ts_to_tsname))
        .join(fpkm, left_on="transcript", right_on="transcript_id(s)", how="left")
    )

    tss_pseudo, pseudo_name = get_pseudogenes(gene, y) if allow_pseudo else (pl.Series(), pl.Series())
    print("Pseudogenes included", pseudo_name)
    pl.DataFrame(dict(transcript=[*tss_all, *tss_pseudo])).write_csv(output / f"{gene}_acceptable_tss.csv")
    ff = y.filter_by_match([*tss_all, *tss_pseudo], match=0.8, match_max=0.8)

    if len(tss_pseudo):
        names_with_pseudo = (
            ff.filter_isin(transcript=tss_pseudo)
            .rename(dict(transcript="maps_to_pseudo"))[["name", "maps_to_pseudo"]]
            .unique("name")
        )
    else:
        names_with_pseudo = pl.DataFrame({"name": ff["name"].unique(), "maps_to_pseudo": ""})
    print(tss_gencode)
    isoforms = (
        ff.filter_isin(transcript=tss_gencode)[["name", "transcript"]]
        .with_columns(isoforms=pl.col("transcript").apply(gtf_all.ts_to_tsname))[["name", "isoforms"]]
        .groupby("name")
        .all()
    )

    ff = GeneFrame(
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
        .filter(~pl.col("seq").apply(lambda x: check_kmers(cast(str, x), kmerset, 18)))
        .filter(~pl.col("seq").apply(lambda x: check_kmers(cast(str, x), trna_rna_kmers, 15)))
        .collect()
        .join(names_with_pseudo, on="name", how="left")
        .join(isoforms, on="name", how="left")
    )
    # print(ff)
    # ff.unique("name").write_parquet(f"output/{gene}_filtered.parquet")
    final = the_filter(ff, overlap=overlap).filter(pl.col("flag") & 16 == 0)
    print(final)
    final.write_parquet(
        output / f"{gene}_final.parquet"
        if overlap < 0
        else output / f"{gene}_final_overlap_{overlap}.parquet"
    )


@click.command()
@click.argument("gene")
@click.option("--output", "-o", type=click.Path(), default="output/")
@click.option("--debug", "-d", is_flag=True)
@click.option("--ignore-revcomp", "-r", is_flag=True)
@click.option("--overlap", "-O", type=int, default=-2)
@click.option("--realign", is_flag=True)
def main(
    gene: str,
    output: str,
    ignore_revcomp: bool,
    overlap: int = -2,
    allow_pseudo: bool = True,
    debug: bool = False,
    realign: bool = False,
):
    if debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        run(
            gene,
            output=output,
            overlap=overlap,
            allow_pseudo=allow_pseudo,
            ignore_revcomp=ignore_revcomp,
            realign=realign,
        )
    except Exception as e:
        raise Exception(f"Failed to run {gene}") from e


if __name__ == "__main__":
    main()
# %%
