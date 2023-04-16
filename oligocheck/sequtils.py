import random
from io import StringIO
from pathlib import Path
from typing import Any

import colorama
import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.axes import Axes

table = str.maketrans("ACGTacgt", "TGCAtgca")


def printc(seq: str):
    for c in seq:
        if c == "A" or c == "a":
            print(colorama.Fore.GREEN + c, end="")
        elif c == "T" or c == "t":
            print(colorama.Fore.RED + c, end="")
        elif c == "C" or c == "c":
            print(colorama.Fore.BLUE + c, end="")
        elif c == "G" or c == "g":
            print(colorama.Fore.YELLOW + c, end="")
        else:
            print(colorama.Fore.WHITE + c, end="")
    print(colorama.Fore.RESET)


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    # https://bioinformatics.stackexchange.com/a/3585
    return seq.translate(table)[::-1]


def pcr(seq: str, primer: str, primer_rc: str) -> str:
    loc = seq.find(primer)
    if loc == -1:
        raise ValueError(f"Primer {primer} not found in sequence {seq}")
    loc_rc = reverse_complement(seq[loc:]).find(primer_rc)
    if loc_rc == -1:
        raise ValueError(f"Primer {primer_rc} not found in sequence {seq}")
    return seq[loc : None if loc_rc == 0 else -loc_rc]


def gen_random_base(n: int) -> str:
    """Generate a random DNA sequence of length n."""
    return "".join(random.choices("ACGT", k=n))


def hamming(seq1: str, seq2: str) -> int:
    """Return the Hamming distance between two sequences."""
    return sum(a != b for a, b in zip(seq1, seq2))


def gc_content(seq: str) -> float:
    return (seq.count("G") + seq.count("C")) / len(seq)


def slide(x: str, n: int = 20):
    return [x[i : i + n] for i in range(len(x) - n + 1)]


def formamide_molar(percent: float = 30) -> float:
    return percent * 10 * 1.13 / 45.04  # density / molar mass


def equal_distance(total: int, choose: int) -> npt.NDArray[np.int_]:
    return np.linspace(0, total - 1, choose).astype(np.int_)


def plot_local_gc_content(ax: Axes, seq: str, window_size: int = 50, **kwargs: Any):
    """Plot windowed GC content on a designated Matplotlib ax."""
    if len(seq) < window_size:
        raise ValueError("Sequence shorter than window size")
    out = np.zeros(len(seq) - window_size + 1, dtype=float)
    curr = gc_content(seq[:window_size]) * window_size
    out[0] = curr
    for i in range(1, len(seq) - window_size):
        curr += (seq[i + window_size - 1] in "GC") - (seq[i - 1] in "GC")
        out[i] = curr
    out /= window_size

    ax.fill_between(np.arange(len(seq) - window_size + 1), out, **(dict(alpha=0.3) | kwargs))  # type: ignore
    ax.set_ylim(bottom=0, top=1)
    ax.set_ylabel("GC (%)")


def parse_sam(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    file_read = [",".join(line.strip().split("\t")[:10]) for line in path.read_text().split("\n")]

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
    return y
