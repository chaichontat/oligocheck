import shlex
import subprocess
from pathlib import Path


def run_bowtie(
    stdin: str | bytes, reference: str, capture_stderr: bool = False, seed_length: int = 17
) -> str:
    return subprocess.run(
        shlex.split(
            # --no-hd No SAM header
            # -k 100 report up to 100 alignments per read
            # -D 20 consecutive seed extension attempts can "fail" before Bowtie 2 moves on
            # -R 3 the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds.
            # -L 17 seed length
            # -i C,2 Seed interval, every 2 bp
            # --score-min G,1,4 f(x) = 1 + 4*ln(read_length)
            f"bowtie2 -x {reference} -U - "
            f"--no-hd -t -k 100 --local -D 20 -R 3 "
            f"-N 0 -L {seed_length} -i C,2 --score-min G,1,4"
        ),
        input=stdin,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE if capture_stderr else None,
        encoding="ascii" if isinstance(stdin, str) else None,
        check=True,
    ).stdout


if __name__ == "__main__":
    out = run_bowtie(Path("zeroth.fastq").read_text(), "data/mm39/mm39", seed_length=17)
