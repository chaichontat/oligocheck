# %%
import json
import logging
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import click
import polars as pl
import primer3
from Levenshtein import distance

from oligocheck.merfish.alignment import run_bowtie
from oligocheck.merfish.external_data import ExternalData, get_rrna
from oligocheck.sequtils import formamide_correction, parse_cigar, parse_sam, slide, tm_hybrid

try:
    profile
except NameError:
    profile = lambda x: x


def tm(s: str) -> float:
    """Approximately the same as the results from the MATLAB script"""
    return primer3.calc_tm(s, mv_conc=300, dv_conc=0, dna_conc=1) + 5


# %%
# GENCODE primary only
gtf = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts.parquet",
    path="data/mm39/gencode.vM32.chr_patch_hapl_scaff.basic.annotation.gtf",
    fasta="data/mm39/combi.fa",
)

gtf_all = ExternalData(
    cache="data/mm39/gencode_vM32_transcripts_all.parquet",
    path="data/mm39/Mus_musculus.GRCm39.109.gtf",
    fasta="data/mm39/combi.fa",
)
# %%
# .set_index("transcript_id")
rrnas = get_rrna("data/mm39/Mus_musculus.GRCm39.ncrna.fa")
trna_rna_kmers = set(
    pl.read_csv("data/mm39/trcombi.txt", separator=" ", has_header=False, new_columns=["kmer", "count"])[
        "kmer"
    ]
)


# %%
@profile
def block_bowtie(gene: str, tss: Iterable[str], temp: Path):
    tss = [f"{gene}_{ts}" for ts in tss]

    for ts in tss:
        if not (temp / f"{ts}.fasta").exists():
            (temp / f"{ts}.fasta").write_text(f">{ts}\n{gtf.get_seq(ts.split('_')[1])}")

    [(temp / f"{ts}.fastq").unlink(missing_ok=True) for ts in tss]
    [(temp / f"{ts}.sam").unlink(missing_ok=True) for ts in tss]

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
                        str(48),
                        "-T",
                        str(55),
                        "-F",
                        "30",
                        "-O",
                    ],
                    check=True,
                ),
                tss,
            )
        )

        # Need some time between python and bowtie2 for files to flush to disk
    while not all((temp / f"{ts}.fastq").exists() for ts in tss):
        print("Waiting for files to flush to disk...")
        time.sleep(0.2)

    with ThreadPoolExecutor() as executor:
        res = executor.map(
            lambda ts: (ts, run_bowtie((temp / f"{ts}.fastq").read_text(), "data/mm39/mm39")),
            tss,
        )

    return list(res)


@profile
def combine_transcripts(sams: Iterable[str]):
    grand = []
    for sam in sams:
        y = parse_sam(sam)
        if len(y.filter(pl.col("transcript") == "*")) > 0:
            raise ValueError("Unmapped reads found")
        grand.append(y)

    # remove duplicates from different bowtie2 runs to the same transcript
    df = (
        pl.concat(grand).unique(["name", "transcript"]).sort(by=["transcript_ori", "pos_start", "transcript"])
    )

    assert (
        df.filter(pl.col("is_ori_seq"))
        .select((pl.col("length") == pl.col("seq").str.n_chars()).alias("is_equal"))["is_equal"]
        .all()
    )
    assert df.null_count().max(axis=1).item() == 0
    return df


@profile
def filter_specifity(
    gene: str,
    tss_all: Iterable[str],
    grand: pl.DataFrame,
    tss_gencode: Iterable[str] | None = None,
    threshold: float = 0.05,
    allow_pseudogene: bool = False,
):
    tsss = set(tss_all)
    tss_gencode = set(tss_gencode) if tss_gencode is not None else set()
    tss_pseudo = set()

    if allow_pseudogene:
        fpkm = pl.read_parquet("data/fpkm/P0_combi.parquet")
        counts = (
            grand.groupby("transcript")
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
            & (
                pl.col("transcript_name").str.starts_with(gene)
                | pl.col("transcript_name").str.starts_with("Gm")
            )
        )

        logging.debug(counts[:50])
        print(f"Including {', '.join(ok['transcript_name'])}.")
        tss_pseudo.update(ok["transcript"])
        assert all(map(lambda x: x.startswith("EN"), tsss))

    # Filter 14-mer rrna/trna
    @profile
    def _filter_specifity(probe: pl.DataFrame):
        if any(x in trna_rna_kmers for x in slide(probe[0, "seq"], n=14)):
            return

        nonspecific = 0
        mapped = set()
        pseudo_binder = False
        for row in probe.iter_rows(named=True):
            transcript = row["transcript"]
            if row["transcript_ori"] == transcript:
                mapped.add(transcript)
                continue  # same transcript as the one designed for

            if transcript in rrnas or "tRNA" in transcript:
                logging.debug(f"Skipping {_} from trna")
                break  # dump this sequence

            if tss_gencode is not None and transcript in tss_gencode:
                # Assuming that it'll bind to said isoform to save computation.
                mapped.add(transcript)
                continue
            elif transcript in tsss:
                # Count only GENCODE transcripts but don't freak out if mapped to a different isoform
                continue
            elif transcript in tss_pseudo:
                pseudo_binder = True
                continue

            if max(parse_cigar(row["cigar"])) > 19:
                # logging.debug(f"Skipping {_} from cigar")
                break

            seq = gtf.get_seq(transcript)

            # To potentially save us some time.
            if distance(seq[row["pos_start"] : row["pos"] + row["pos_end"] + 2], row["seq"]) < 5:
                # logging.debug(f"Skipping {_} from distance")
                break

            # ns = nonspecific_test(
            #     row["seq"],
            #     seq[max(0, row["pos"] - 2) : min(row["pos"] + len(row["seq"]) + 2, len(["seq"]))],
            #     47,
            # )[1]["bound"]

            # if not ontarget_any and ns < threshold:
            #     nonspecific = max(ns, nonspecific)

            # if not ontarget_any and ns >= threshold:
            #     logging.debug(f"Skipping {_} from offtarget")
            #     break

        else:
            if tss_gencode is not None:
                assert len(mapped) > 0
                probe = probe.with_columns(pl.lit(len(mapped)).alias("n_mapped"))
            probe = probe.with_columns(
                nonspecific_binding=pl.lit(nonspecific), pseudobinder=pl.lit(pseudo_binder)
            )
            # logging.debug(f"Keeping {len(probe)}")
            return probe

    res = []
    for _, rows in grand.groupby("name"):
        df = _filter_specifity(rows)
        if df is not None:
            res.append(df)

    if not res:
        print(f"No candidates found for {gene}")
        return pl.DataFrame()
    return pl.concat(res), list(tsss)


# def filter_candidates(res: pd.DataFrame, stringency: Stringency, n_cand: int = 300):
#     if stringency.overlap > 0 or len(res) < n_cand:
#         # Not dropping anything here since we don't have enough.
#         return res

#     curr = res.n_mapped.max()
#     picked = pd.DataFrame()
#     # shuffle to get uniform coverage in case we get to n_cand early.
#     while len(picked) < n_cand:
#         shuffled = res[res.n_mapped == curr].sample(frac=1)
#         if len(shuffled) + len(picked) > n_cand:
#             shuffled = shuffled.iloc[: n_cand - len(picked)]
#         picked = pd.concat([picked, shuffled]) if len(picked) else shuffled
#     return picked


@profile
def calc_thermo(picked: pl.DataFrame):
    return picked.with_columns(
        [
            pl.col("seq")
            .apply(
                lambda x: primer3.calc_hairpin_tm(x, mv_conc=390, dv_conc=0, dntp_conc=0, dna_conc=1)
                + formamide_correction(x)
            )
            .alias("hp"),
            pl.col("seq")
            .apply(
                lambda x: primer3.calc_homodimer_tm(x, mv_conc=390, dv_conc=0, dntp_conc=0, dna_conc=1)
                + formamide_correction(x)
            )
            .alias("homodimer"),
            pl.col("seq").apply(tm_hybrid).alias("tm"),
        ]
    )


# %%
# gene = "Pclaf"


@profile
def run(gene: str):
    tss_gencode = gtf.filter_gene(gene)["transcript_id"]
    longest = max(tss_gencode, key=lambda x: len(gtf.get_seq(x)))
    tss_all = gtf_all.filter_gene(gene)["transcript_id"]

    print(f"Running {gene} with {len(tss_gencode)} transcripts")

    temp = Path("temp")
    temp.mkdir(exist_ok=True)

    sams = block_bowtie(gene, [longest], temp)
    grand = combine_transcripts(map(lambda x: x[1], sams))
    print(f"{gene}: Found {len(grand)} SAM entries")
    res, allowed_transcripts = filter_specifity(
        gene, tss_all, grand, tss_gencode=tss_gencode, allow_pseudogene=True
    )

    res = (
        res.lazy()
        .with_columns(
            [
                pl.col("transcript").apply(gtf_all.eid_to_ts).alias("transcript_name"),
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
            ]
        )
        .collect()
    )

    if not len(res):
        raise ValueError(f"No candidates found for {gene}")
    print(f"{gene}: Found {len(res.unique('name'))} probes after specificity filtering")

    # %%
    res = res.sort(
        by=["is_ori_seq", "transcript_ori", "pos_start", "transcript"], descending=[True, False, False, False]
    ).unique(subset=["name"], keep="first")

    picked = calc_thermo(res)
    return picked.drop(columns=["rnext", "pnext", "tlen"]), dict(
        gene=gene,
        generated_on=tss_gencode.to_list(),
        allowed=allowed_transcripts,
        allowed_name=list(map(gtf_all.ts_to_gene, allowed_transcripts)),
    )


# %%


@click.command()
@click.argument("gene")
@click.option("--output", "-o", type=click.Path(), default="output/")
@click.option("--debug", "-d", is_flag=True)
def main(gene: str, output: str | Path = "output/", debug: bool = False):
    if debug:
        logging.basicConfig(level=logging.DEBUG)
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        m, md = run(gene)
    except Exception as e:
        raise Exception(f"Failed to run {gene}") from e

    Path(output).mkdir(exist_ok=True, parents=True)
    m.write_parquet(f"{output}/{gene}.parquet")
    with open(f"{output}/{gene}.json", "w") as f:
        json.dump(md, f)


if __name__ == "__main__":
    main()


# %%

# -L ignore kmer < 2 c is counter bit. Want to be such that most kmers use only 1 counter, so 8.
# jellyfish count -m 18 -o output -t 32 -s 10G -L 2 -c 3 combi.fa
# jellyfish dump -c -L 2 output > output.txt
# count = dict()
# for i, x in jsorted.iterrows():
#     count[x[0]] = x[1]

# %%
# %%
