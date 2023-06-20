# %%
import shlex
import subprocess
from pathlib import Path
from typing import Collection

from oligocheck.boilerplate import check_if_posix
from oligocheck.merfish.alignment import gen_fasta


@check_if_posix
def run_jellyfish(
    seqs: Collection[str] | str,
    out: str | Path,
    kmer: int,
    *,
    hash_size: str = "10G",
    minimum: int = 1,
    counter: int = 2,
    thread: int = 16,
    both_strands: bool = False,
):
    """
    http://www.cs.cmu.edu/~ckingsf/software/jellyfish/jellyfish-manual-1.1.pdf
    """
    if isinstance(seqs, str):
        stdin = Path(seqs).read_text()
    else:
        stdin = gen_fasta(map(str, range(len(seqs))), seqs).getvalue()

    res = subprocess.run(
        shlex.split(
            rf"jellyfish count  -o /dev/stdout -m {kmer} -t {thread} -s {hash_size} -L {minimum} -c {counter} {'--both-strands' if both_strands else ''} /dev/stdin"
        ),
        input=stdin,
        stdout=subprocess.PIPE,
        encoding="ascii",
        check=True,
    )
    res = subprocess.run(
        # -c == column format
        shlex.split(rf"jellyfish dump -c -L 1 /dev/stdin > {out}"),
        input=res.stdout,
        stdout=subprocess.PIPE,
        encoding="ascii",
        check=True,
    )
