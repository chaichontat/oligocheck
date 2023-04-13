# %%
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Iterable, Literal

import click
import pandas as pd
import primer3

from oligocheck.merfish.external_data import all_transcripts, gene_to_eid, get_gencode, get_rrna, get_seq
from oligocheck.merfish.nnupack import nonspecific_test
from oligocheck.sequtils import formamide_molar, slide


def tm(s: str) -> float:
    """Approximately the same as the results from the MATLAB script"""
    return primer3.calc_tm(s, mv_conc=300, dv_conc=0, dna_conc=1) + 5


@dataclass(frozen=True)
class Stringency:
    min_tm: float = 47
    max_tm: float = 52
    min_gc: float = 35
    max_gc: float = 65
    overlap: float = -1
    hairpin_tm: float = 30
    unique: bool = False  # to the extent that bowtie2 cannot detect


# fmt: off
@dataclass(frozen=True)
class Stringencies:
    # high   = Stringency(min_tm=49, min_gc=35, max_gc=65, hairpin_tm=40, filter_quad_c=True)
    high   = Stringency(min_tm=49, min_gc=35, max_gc=65, hairpin_tm=35, unique=False)
    medium = Stringency(min_tm=49, min_gc=25, max_gc=70, hairpin_tm=35, unique=False)
    low    = Stringency(min_tm=48, min_gc=25, max_gc=70, hairpin_tm=35, unique=False)
    high_ol   = Stringency(min_tm=48, min_gc=35, max_gc=65, hairpin_tm=20, unique=False, overlap=15)
    medium_ol = Stringency(min_tm=48, min_gc=25, max_gc=70, hairpin_tm=30, unique=False, overlap=15)
    low_ol    = Stringency(min_tm=48, min_gc=25, max_gc=70, hairpin_tm=35, unique=False, overlap=15)
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
                        "blockParse.py",
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
                f"bowtie2 -x ../mm39/mm39 {(temp / f'{ts}.fastq').as_posix()} -t --no-hd -k 2 --local -D 20 -R 3 -N 1 -L 20 -i C,4 --score-min G,1,4 -S {(temp / f'{ts}.sam').as_posix()}",
                shell=True,
                check=True,
            ),
            tss,
        )

    while not all((temp / f"{ts}.sam").exists() for ts in tss):
        print("Waiting for files to flush to disk...")
        time.sleep(0.2)


def combine_transcripts(gene: str, tss: Iterable[str], temp: Path):
    grand = []
    for ts in [f"{gene}_{ts}" for ts in tss]:
        file_read = [
            ",".join(line.strip().split("\t")[:10]) for line in (temp / f"{ts}.sam").read_text().split("\n")
        ]

        y = pd.read_csv(
            StringIO("\n".join(file_read)),
            sep=",",
            header=None,
            names=[
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
            ],
        )
        y.transcript = y.transcript.apply(lambda x: x.split(".")[0])
        if "*" in y.transcript.value_counts().index:
            raise ValueError("Unmapped reads found")
        grand.append(y)

    grand = pd.concat(grand).drop_duplicates(subset=["transcript", "seq"], keep="first")
    return grand


def filter_specifity(
    tss_all: Iterable[str],
    tss_gencode: Iterable[str],
    grand: pd.DataFrame,
    stringency: Stringency,
    threshold: float = 0.05,
):
    grand["ok"] = False
    grand["nonspecific_binding"] = 0.0
    tsss = set(tss_all)
    tssg = set(tss_gencode)

    grand["ok"] = False
    grand.loc[grand.transcript.isin(tsss), "ok"] = True
    res = []

    for _, rows in grand.groupby("seq"):
        # Filter 14-mer rrna/trna
        if any(x in trna_rna_kmers for x in slide(rows.iloc[0].seq, n=14)):
            continue

        all_ontarget = rows.transcript.isin(tsss).all()
        if stringency.unique and not all_ontarget:
            continue

        nonspecific = 0
        for _, row in rows.iterrows():
            if all_ontarget:
                continue  # No need to test. Cannot break because of else.

            if row.transcript in rrnas or "tRNA" in row.transcript:
                break

            # Actual nonspecific test
            seq = get_seq(row.transcript)
            ns = nonspecific_test(
                row.seq,
                seq[max(0, row.pos - 3) : min(row.pos + len(row.seq) + 3, len(seq))],
                37,
            )[1]["bound"]
            if ns > threshold:
                break
            nonspecific = max(ns, nonspecific)
        else:
            # Counting only GENCODE transcripts
            rows["n_mapped"] = len(rows[rows.transcript.isin(tssg)])
            rows["nonspecific_binding"] = nonspecific
            res.append(rows)

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
            lambda x: primer3.calc_hairpin_tm(x, mv_conc=300, dv_conc=0, dntp_conc=0, dna_conc=1)
            - 0.65 * 30  # formamide
        ),
        homodimer=picked["seq"].map(
            lambda x: primer3.calc_homodimer_tm(x, mv_conc=300, dv_conc=0, dntp_conc=0, dna_conc=1)
            - 0.65 * 30  # formamide
        ),
        tm=picked["seq"].map(
            lambda x: primer3.calc_tm(
                x,
                mv_conc=300,
                dv_conc=0,
                dna_conc=1,
                formamide_conc=formamide_molar(0.3),
            )
        ),
    )


def filter_thermo(picked: pd.DataFrame, stringency: Stringency):
    return picked[(picked.hp < stringency.hairpin_tm) & (picked.homodimer < 40)]


# %%
# stg = "low"
# gene = "B2m"


def run(gene: str, stg: Literal["high", "medium", "low", "low_ol"]):
    stringency: Stringency = getattr(Stringencies, stg)

    eid = gene_to_eid(gene)
    tss_gencode = tuple(everything[everything.gene_id == eid]["transcript_id"].values)
    tss_all = all_transcripts(gene)

    print(f"Running {gene} with {len(tss_gencode)} transcripts")

    temp = Path("temp")
    temp.mkdir(exist_ok=True)

    block_bowtie(gene, tss_gencode, stringency, temp)
    grand = combine_transcripts(gene, tss_gencode, temp)
    res = filter_specifity(tss_all, tss_gencode, grand, stringency)
    res = res.drop_duplicates(subset=["seq"], keep="first")
    # picked = filter_candidates(res, stringency)
    picked = calc_thermo(res)
    picked = filter_thermo(picked, stringency)

    return picked.sort_values(by=["transcript", "pos"]).drop(
        columns=["mapq", "cigar", "flag", "rnext", "pnext", "tlen", "ok"]
    )


@click.command()
@click.argument("gene")
@click.argument("stringency")
@click.option("--output", "-o", type=click.Path())
def main(gene: str, stringency: str, output: str):
    m = run(gene, stringency)
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
# %%
# %%
# %%
# %%
# %%
# %%
# %%
# %%
