import io
import shlex
import subprocess
from pathlib import Path

import polars as pl
from Bio import AlignIO


def gen_fastq(names: pl.Series, seqs: pl.Series) -> io.StringIO:
    f = io.StringIO()
    for name, seq in zip(names, seqs):
        f.write(f"@{name}\n")
        f.write(seq + "\n")
        f.write("+\n")
        f.write("~" * len(seq) + "\n")
    return f


def gen_fasta(names: pl.Series, seqs: pl.Series) -> io.StringIO:
    f = io.StringIO()
    for name, seq in zip(names, seqs):
        f.write(f">{name}\n")
        f.write(seq + "\n")
    return f


def gen_bowtie_index(stdin: str, path: str, name: str) -> bytes:
    (Path(path) / (name + ".fasta")).write_text(stdin)
    return subprocess.run(
        shlex.split(f"bowtie2-build {path}/{name}.fasta {path}/{name}"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    ).stdout


def run_bowtie(
    stdin: str | bytes,
    reference: str,
    capture_stderr: bool = False,
    seed_length: int = 15,
    n_return: int = 100,
    threads: int = 16,
    threshold: int = 15,
) -> str:
    return subprocess.run(
        shlex.split(
            # A base that matches receives a bonus of +2 be default.
            # A mismatched base at a high-quality position in the read receives a penalty of -6 by default.
            # --no-hd No SAM header
            # -k 100 report up to 100 alignments per read
            # -D 20 consecutive seed extension attempts can "fail" before Bowtie 2 moves on
            # -R 3 the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds.
            # -L 17 seed length
            # -i C,2 Seed interval, every 2 bp
            # --score-min G,1,4 f(x) = 1 + 4*ln(read_length)
            # --score-min L,0,-0.6 f(x) = -0.6*read_length
            f"bowtie2 -x {reference} -U - "
            f"--no-hd -t {f'-k {n_return}' if n_return > 0 else '-a'} --local -D 20 -R 3 "
            f"--score-min L,{threshold*2},0 --mp 1,1 --ignore-quals "
            f"-N 0 -L {seed_length} -i C,2 -p {threads}"
        ),
        input=stdin,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE if capture_stderr else None,
        encoding="ascii" if isinstance(stdin, str) else None,
        check=True,
    ).stdout


def run_mafft(stdin: str | bytes) -> str:
    return subprocess.run(
        ["mafft-linsi", "--op", "2", "-"],
        input=stdin,
        encoding="ascii" if isinstance(stdin, str) else None,
        capture_output=True,
        check=True,
    ).stdout


def parse_mafft(s: str) -> AlignIO.MultipleSeqAlignment:
    return AlignIO.read(io.StringIO(s), "fasta")


if __name__ == "__main__":
    out = run_bowtie(Path("zeroth.fastq").read_text(), "data/mm39/mm39", seed_length=17)
