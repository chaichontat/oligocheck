# %%
import logging
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Iterable

import click
import polars as pl

from oligocheck.algorithms import find_overlap, find_overlap_weighted
from oligocheck.merfish.alignment import run_bowtie
from oligocheck.merfish.external_data import ExternalData, get_rrna
from oligocheck.seqcalc import hp_fish, tm_fish

# from oligocheck.merfish.filtration import the_filter
from oligocheck.sequtils import parse_sam

sys.setrecursionlimit(1500)
pl.Config.set_fmt_str_lengths(100)
pl.Config.set_tbl_rows(30)

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
    tss = [f"{gene}_{ts}" for ts in tss]

    for ts in tss:
        if not (temp / f"{ts}.fasta").exists():
            (temp / f"{ts}.fasta").write_text(f">{ts}\n{gtf.get_seq(ts.split('_')[1])}")

    [(temp / f"{ts}.fastq").unlink(missing_ok=True) for ts in tss]

    with ThreadPoolExecutor() as executor:
        list(
            executor.map(
                lambda ts: subprocess.run(
                    [
                        "python",
                        "oligocheck/merfish/blockParse.py",
                        "-f",
                        (temp / f"{ts}.fasta").as_posix(),
                        "-t",
                        str(50),
                        "-T",
                        str(52),
                        "-F",
                        "30",
                        "-O",
                    ],
                    check=True,
                ),
                tss,
            )
        )
    while not all((temp / f"{ts}.fastq").exists() for ts in tss):
        print("Waiting for files to flush to disk...")
        time.sleep(0.2)

    with ThreadPoolExecutor() as executor:
        res = executor.map(
            lambda ts: (
                ts,
                run_bowtie(
                    (temp / f"{ts}.fastq").read_text(),
                    "data/mm39/mm39",
                    seed_length=12,
                    threshold=16,
                    n_return=-1,
                ),
            ),
            tss,
        )

    return list(res)


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


def filter_match(y: pl.DataFrame, acceptable_tss: Iterable[str], match: int = 21, match_max: int = 16):
    nogo = y.filter(
        (pl.col("match").gt(match))
        & (pl.col("match_max").gt(match_max))
        & ~pl.col("transcript").str.contains("|".join([f"({x})" for x in acceptable_tss]))
    )
    return y.filter(~pl.col("name").is_in(nogo["name"]))


# %%


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
    df = ddf.filter(pl.col("priority") > 0).collect()

    selected_global = set()
    for ts in tss:
        this_transcript = df.filter(pl.col("transcript_ori") == ts)
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
            if i == 1:
                ols = find_overlap(run["pos_start"], run["pos_end"], overlap=overlap)
            else:
                ols = find_overlap_weighted(run["pos_start"], run["pos_end"], priorities, overlap=overlap)
            sel_local = set(run[ols]["index"].to_list())
            print(i, len(sel_local))
            if len(sel_local) > n:
                break

        selected_global |= sel_local  # type: ignore

    return df.filter(pl.col("index").is_in(selected_global))


# %%
# %%

# %%


def the_filter(df: pl.DataFrame, overlap: int = -1):
    out = []
    for _, group in df.groupby("gene"):
        out.append(
            handle_overlap(
                gtf,
                group,
                criteria=[
                    # fmt: off
                    pl.col("tm").is_between(50, 54) & (pl.col("oks") > 4) & (pl.col("hp") < 35),
                    pl.col("tm").is_between(50, 54) & (pl.col("oks") > 3) & (pl.col("hp") < 35),
                    pl.col("tm").is_between(49, 56) & (pl.col("oks") > 3) & (pl.col("hp") < 40),
                    pl.col("tm").is_between(49, 56) & (pl.col("oks") > 2) & (pl.col("hp") < 40),
                    pl.col("tm").is_between(49, 56) & (pl.col("oks") > 2) & (pl.col("hp") < 45),
                    # fmt: on
                ],
                overlap=overlap,
            )
        )
    return pl.concat(out)


# %%
def run(gene: str, overlap: int = -2, allow_pseudo: bool = False):
    Path("temp").mkdir(exist_ok=True)
    # wants = list(filter(lambda x: x, Path("celltypegenes.csv").read_text().splitlines()))
    # gene = wants[7]
    tss_gencode = gtf.filter_gene(gene)["transcript_id"]
    tss_all = gtf_all.filter_gene(gene)["transcript_id"]
    longest = max(tss_gencode, key=lambda x: len(gtf.get_seq(x)))

    y = parse_sam(block_bowtie(gene, [longest], Path("temp"))[0][1])
    y.write_parquet(f"output/{gene}_all.parquet")
    tss_pseudo = get_pseudogenes(gene, y) if allow_pseudo else pl.Series()
    offtargets = (
        y["transcript"]
        .value_counts()
        .sort("counts", descending=True)
        .with_columns(name=pl.col("transcript").apply(gtf_all.ts_to_gene))
    )

    offtargets.write_parquet(f"output/{gene}_offtargets.parquet")

    y = y.join(
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
    )

    ff = (
        filter_match(y, [*tss_all, *tss_pseudo], match=21, match_max=17)
        .lazy()
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
                pl.col("seq").apply(tm_fish).alias("tm") - 0.65 * 30,
                pl.col("seq").apply(hp_fish).alias("hp") - 0.65 * 30,
            ]
        )
        .with_columns(oks=pl.sum(pl.col("^ok_.*$")))
        .collect()
    )

    final = the_filter(ff, overlap=overlap)
    final.write_parquet(f"output/{gene}_final.parquet")


@click.command()
@click.argument("gene")
# @click.option("--output", "-o", type=click.Path(), default="output/")
@click.option("--debug", "-d", is_flag=True)
@click.option("--allow-pseudo", "-p", is_flag=True)
@click.option("--overlap", "-O", type=int, default=-2)
def main(gene: str, overlap: int = -2, allow_pseudo: bool = False, debug: bool = False):
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
    main()
