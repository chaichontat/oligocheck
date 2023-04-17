# %%
import logging
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal

import click
import pandas as pd
import primer3
from Levenshtein import distance

from oligocheck.merfish.external_data import (
    all_transcripts,
    gene_to_eid,
    get_gencode,
    get_rrna,
    get_seq,
)
from oligocheck.merfish.nnupack import nonspecific_test
from oligocheck.sequtils import (
    formamide_correction,
    gc_content,
    parse_sam,
    reverse_complement,
    slide,
    tm_hybrid,
)

# logging.basicConfig(level=logging.DEBUG)
# logging.getLogger().setLevel(logging.DEBUG)


def tm(s: str) -> float:
    """Approximately the same as the results from the MATLAB script"""
    return primer3.calc_tm(s, mv_conc=300, dv_conc=0, dna_conc=1) + 5


@dataclass(frozen=True)
class Stringency:
    min_tm: float = 47
    max_tm: float = 57
    min_gc: float = 35
    max_gc: float = 65
    overlap: float = -1
    hairpin_tm: float = 30
    unique: bool = False  # to the extent that bowtie2 cannot detect
    max_notok: int = 2


# fmt: off
@dataclass(frozen=True)
class Stringencies:
    # high   = Stringency(min_tm=49, min_gc=35, max_gc=65, hairpin_tm=40, filter_quad_c=True)
    high   = Stringency(min_tm=49, min_gc=35, max_gc=65, hairpin_tm=30, unique=False)
    medium = Stringency(min_tm=49, min_gc=25, max_gc=70, hairpin_tm=45, unique=False)
    low    = Stringency(min_tm=48, min_gc=25, max_gc=70, hairpin_tm=35, unique=False)
# fmt: on


# %%
# GENCODE primary only
everything = get_gencode("data/mm39/gencode_vM32_transcripts.parquet")
rrnas = get_rrna("data/mm39/rrna.fa")
trna_rna_kmers = set(
    pd.read_csv("data/mm39/trcombi.txt", sep=" ", header=None, names=["counts"], index_col=0)["counts"].index
)


# %%
# gene = 'Pclaf'
# stg = 'high'
def block_bowtie(gene: str, tss: Iterable[str], stringency: Stringency, temp: Path):
    tss = [f"{gene}_{ts}" for ts in tss]

    for ts in tss:
        if not (temp / f"{ts}.fasta").exists():
            (temp / f"{ts}.fasta").write_text(f">{ts}\n{get_seq(ts.split('_')[1])}")

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
                        str(stringency.min_tm),
                        "-T",
                        str(stringency.max_tm),
                        "-F",
                        "30",
                        "-O",
                    ],
                    # + (
                    #     ["-S", str(-stringency.overlap)]
                    #     if stringency.overlap <= 0
                    #     else []
                    # )
                    # + (["-O"] if stringency.overlap > 0 else []),
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
        executor.map(
            lambda ts: subprocess.run(
                # --no-hd No SAM header
                # -a report all alignments
                # -D 20 consecutive seed extension attempts can "fail" before Bowtie 2 moves on
                # -R 3 the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds.
                # -L 17 seed length
                # -i C,2 Seed interval, every 2 bp
                # --score-min G,1,4 f(x) = 1 + 4*ln(read_length)
                f"bowtie2 -x data/mm39/mm39 {(temp / f'{ts}.fastq').as_posix()} "
                f"--no-hd -a --local -D 20 -R 3 -N 0 -L 17 -i C,2 --score-min G,1,4 "
                f"-S {(temp / f'{ts}.sam').as_posix()}",
                shell=True,
                check=True,
            ),
            tss,
        )

    while not all((temp / f"{ts}.sam").exists() for ts in tss):
        print("Waiting for files to flush to disk...")
        time.sleep(0.2)


name_splitter = re.compile(r"(.+)_(.+):(\d+)-(\d+)")


def filter_oks(df: pd.DataFrame, max_notok: int = 2):
    oks = df.columns[df.columns.str.contains("ok_")]
    df["oks"] = df[oks].sum(axis=1)
    return df[df.oks > len(oks) - max_notok]


def cSpecStackFilter(seq: str):
    seq = reverse_complement(seq)
    for i in range(6):
        if seq[i : i + 6].count("C") >= 4:
            return False
    return True


def ACompfilter(seq: str):
    return seq.count("T") / len(seq) < 0.28


def combine_transcripts(gene: str, tss: Iterable[str], temp: Path):
    grand = []
    for ts in [f"{gene}_{ts}" for ts in tss]:
        y = parse_sam(temp / f"{ts}.sam")
        if "*" in y.transcript.value_counts().index:
            raise ValueError("Unmapped reads found")
        grand.append(y)

    # remove duplicates from different bowtie2 runs to the same transcript
    grand = pd.concat(grand).drop_duplicates(subset=["name", "transcript"])
    df = grand.sort_values(by=["transcript_ori", "pos_start", "transcript"])

    assert (df[df.is_ori_seq].length == df[df.is_ori_seq].seq.map(len)).all()
    assert not df.isna().any().any()
    return df


def filter_specifity(
    tss_all: Iterable[str],
    grand: pd.DataFrame,
    tss_gencode: Iterable[str] | None = None,
    threshold: float = 0.05,
):
    grand["nonspecific_binding"] = 0.0
    tsss = set(tss_all)
    tss_gencode = set(tss_gencode) if tss_gencode is not None else set()

    # Filter 14-mer rrna/trna
    def _filter_specifity(rows: pd.DataFrame):
        if any(x in trna_rna_kmers for x in slide(rows.iloc[0].seq, n=14)):
            return

        # if stringency.unique and not all_ontarget:
        #     continue
        nonspecific = 0
        mapped = set()
        for _, row in rows.iterrows():
            if row["transcript_ori"] == row.transcript:
                mapped.add(row.transcript)
                continue  # same transcript as the one designed for

            if row.transcript in rrnas or "tRNA" in row.transcript:
                logging.debug(f"Skipping {_} from trna")
                break  # dump this sequence

            # Count only GENCODE transcripts but don't freak out if mapped to a diffrent isoform
            ontarget_any = row.transcript in tsss

            # Assuming that it'll bind to said isoform to save computation.
            if tss_gencode is not None and row.transcript in tss_gencode:
                mapped.add(row.transcript)
                continue
            elif ontarget_any:
                continue

            seq = get_seq(row.transcript)

            # To potentially save us some time.
            if distance(seq[row.pos : row.pos + len(row.seq)], row.seq) < 5:
                logging.debug(f"Skipping {_} from distance")
                break

            ns = nonspecific_test(
                row.seq,
                seq[max(0, row.pos - 3) : min(row.pos + len(row.seq) + 3, len(seq))],
                47,
            )[1]["bound"]

            if not ontarget_any and ns < threshold:
                nonspecific = max(ns, nonspecific)

            if not ontarget_any and ns >= threshold:
                logging.debug(f"Skipping {_} from offtarget")
                break

        else:
            if tss_gencode is not None:
                assert len(mapped) > 0
                rows["n_mapped"] = len(mapped)
            rows["nonspecific_binding"] = nonspecific
            logging.debug(f"Keeping {len(rows)}")
            return rows

    res = []
    for _, rows in grand.groupby("name"):
        _ = _filter_specifity(rows)
        if _ is not None:
            res.append(_)

    # with ThreadPoolExecutor(6) as executor:
    #     tosplit = list(grand.groupby("name"))
    #     # split tosplit uniformly into 6 chunks
    #     chunks = [tosplit[i::6] for i in range(6)]

    #     futures.append(executor.submit(_filter_specifity, rows))
    #     for future in as_completed(futures):
    #         _ = future.result()
    #         if _ is not None:
    #             res.append(_)
    return pd.concat(res)


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


def calc_thermo(picked: pd.DataFrame):
    return picked.assign(
        hp=picked["seq"].map(
            lambda x: primer3.calc_hairpin_tm(x, mv_conc=390, dv_conc=0, dntp_conc=0, dna_conc=1)
            + formamide_correction(x)
        ),
        homodimer=picked["seq"].map(
            lambda x: primer3.calc_homodimer_tm(x, mv_conc=390, dv_conc=0, dntp_conc=0, dna_conc=1)
            + formamide_correction(x)
        ),
        tm=picked["seq"].map(
            tm_hybrid
            # primer3.calc_tm(
            #     x,
            #     mv_conc=300,
            #     dv_conc=0,
            #     dna_conc=1,
            #     formamide_conc=formamide_molar(0.3),
            # )
        ),
    )


def filter_thermo(picked: pd.DataFrame, stringency: Stringency):
    return picked[(picked.hp < stringency.hairpin_tm) & (picked.homodimer < 40)]


# %%
stg = "high"
gene = "Pclaf"


def run(gene: str, stg: Literal["high", "medium", "low", "low_ol"]):
    stringency: Stringency = getattr(Stringencies, stg)

    eid = gene_to_eid(gene)
    tss_gencode = everything[(everything.gene_id == eid) & everything.transcript_support_level == 1][
        "transcript_id"
    ]
    tss_all = all_transcripts(gene)

    print(f"Running {gene} with {len(tss_gencode)} transcripts")

    temp = Path("temp")
    temp.mkdir(exist_ok=True)

    block_bowtie(gene, tss_gencode, stringency, temp)
    grand = combine_transcripts(gene, tss_gencode, temp)
    print(f"{gene}: Found {len(grand)} SAM entries")

    grand["ok_quad_c"] = ~grand["seq"].str.contains("GGGG")
    grand["ok_quad_a"] = ~grand["seq"].str.contains("TTTT")
    grand["ok_comp_a"] = grand["seq"].map(ACompfilter)
    grand["ok_stack_c"] = grand["seq"].map(cSpecStackFilter)
    grand["gc_content"] = grand["seq"].map(gc_content)
    grand["ok_gc"] = (grand.gc_content >= 0.35) & (grand.gc_content <= 0.65)
    # grand = filter_oks(grand, stringency.max_notok)
    print(f"{gene}: Found {len(grand)} SAM entries after filtering")

    res = filter_specifity(tss_all, grand, tss_gencode=tss_gencode)
    res = res.sort_values(
        by=["is_ori_seq", "transcript_ori", "pos_start", "transcript"], ascending=[False, True, True, True]
    ).drop_duplicates(subset=["name"], keep="first")

    # picked = filter_candidates(res, stringency)
    picked = calc_thermo(res)
    # picked = filter_thermo(picked, stringency)
    return picked.drop(columns=["mapq", "cigar", "flag", "rnext", "pnext", "tlen"])


# %%


@click.command()
@click.argument("gene")
@click.argument("stringency")
@click.option("--output", "-o", type=click.Path())
def main(gene: str, stringency: str, output: str):
    try:
        m = run(gene, stringency)
    except Exception as e:
        raise Exception(f"Failed to run {gene} with {stringency}") from e
    Path(output).mkdir(exist_ok=True, parents=True)
    m.to_parquet(f"{output}/{gene}_{stringency}.parquet", index=False)


if __name__ == "__main__":
    main()


# %%


# %%


# return selected, filtered

# %%
# if __name__ == "__main__":
# selected, filtered = run("Pclaf", "high")


# %%


# picked = filter_candidates(res, stringency)
# picked = calc_thermo(picked, stringency)
# filtered = picked[(picked.hp < stringency.hairpin_tm) & (picked.homodimer < 55)]
# selected = handle_overlap(tss, filtered, stringency)

# if stringency.filter_quad_c:
#     selected = selected[~selected['seq'].str.contains("CCCC")]


# print(len(selected))
# %%

# -L ignore kmer < 2 c is counter bit. Want to be such that most kmers use only 1 counter, so 8.
# jellyfish count -m 18 -o output -t 32 -s 10G -L 2 -c 3 combi.fa
# jellyfish dump -c -L 2 output > output.txt
# count = dict()
# for i, x in jsorted.iterrows():
#     count[x[0]] = x[1]

# %%


# %%


# %%


# %%


# %%
# %%


# %%


# %%
