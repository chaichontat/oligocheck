# %%
import random
import re
from io import StringIO
from typing import Any

import colorama
import numpy as np
import numpy.typing as npt
import polars as pl
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from matplotlib.axes import Axes

from oligocheck.merfish.pairwise import pairwise_alignment

table = str.maketrans("ACGTacgt ", "TGCAtgca ")
name_splitter = re.compile(r"(.+)_(.+):(\d+)-(\d+)")

c = re.compile(r"(\d+)S(\d+)M")
c2 = re.compile(r"(\d+)M")


def parse_cigar(s: str, m_only: bool = False) -> tuple[int, ...]:
    if not m_only:
        try:
            return tuple(map(int, c.findall(s)[0]))
        except IndexError:
            try:
                return 0, max(map(int, (c2.findall(s))))
            except ValueError:
                return 0, 0
            except IndexError:
                return 0, 0
    return tuple(map(int, c2.findall(s)))


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


def equal_distance(total: int, choose: int) -> npt.NDArray[np.int_]:
    return np.linspace(0, total - 1, choose).astype(np.int_)


def plot_gc_content(ax: Axes, seq: str, window_size: int = 50, **kwargs: Any):
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


def parse_sam(sam: str, split_name: bool = True) -> pl.DataFrame:
    # file_read = [",".join(line.strip().split("\t")[:10]) for line in sam.split("\n")]

    # s = (
    #     pl.DataFrame(dict(strs=[sam]))
    #     .lazy()
    #     .with_columns(pl.col("strs").str.split("\n"))
    #     .explode("strs")
    #     .with_columns(pl.col("strs").str.strip().str.split("\t").arr.slice(0, 10))
    #     .with_row_count("id")
    #     .explode("strs")
    #     .with_columns(col_nm="string_" + pl.arange(0, pl.count()).cast(pl.Utf8).str.zfill(2).over("id"))
    #     .sort(["col_nm", "id"])
    #     .collect()
    # )

    # faster than pivot
    key_optional = {
        "AS": "aln_score",
        "XS": "aln_score_best",
        "XN": "n_ambiguous",
        "XM": "n_mismatches",
        "XO": "n_opens",
        "XG": "n_extensions",
        "NM": "edit_distance",
        "MD": "mismatched_reference",
        "YT": "pair_state",
    }

    df = (
        pl.read_csv(StringIO(sam), separator="\n", has_header=False)
        .lazy()
        .with_row_count("id")
        .with_columns(temp=pl.col("column_1").str.split_exact("\t", 9))
        .unnest("temp")
        .rename(
            {
                f"field_{i}": x
                for i, x in enumerate(
                    ["name", "flag", "transcript", "pos", "mapq", "cigar", "rnext", "pnext", "tlen", "seq"]
                )
            }
        )
        .with_columns(
            flag=pl.col("flag").cast(pl.UInt16),
            pos=pl.col("pos").cast(pl.UInt32),
            mapq=pl.col("mapq").cast(pl.UInt8),
            aln_score=pl.col("column_1").str.extract(r"AS:i:(\d+)").cast(pl.UInt16),
            aln_score_best=pl.col("column_1").str.extract(r"XS:i:(\d+)").cast(pl.UInt16),
            n_ambiguous=pl.col("column_1").str.extract(r"XN:i:(\d+)").cast(pl.UInt16),
            n_mismatches=pl.col("column_1").str.extract(r"XM:i:(\d+)").cast(pl.UInt16),
            n_opens=pl.col("column_1").str.extract(r"XO:i:(\d+)").cast(pl.UInt16),
            n_extensions=pl.col("column_1").str.extract(r"XG:i:(\d+)").cast(pl.UInt16),
            edit_distance=pl.col("column_1").str.extract(r"NM:i:(\d+)").cast(pl.UInt16),
            mismatched_reference=pl.col("column_1").str.extract(r"MD:Z:(\S+)"),
        )
        .drop(["column_1", "mapq", "rnext", "pnext", "tlen"])
        .with_columns(
            [
                pl.when(pl.col("transcript").str.contains(r"(.*)\.\d+"))
                .then(pl.col("transcript").str.extract(r"(.*)\.\d+"))
                .otherwise(pl.col("transcript"))
                .alias("transcript")
            ]
            + [
                pl.col("name").str.extract(r"(.+)_(.+):(\d+)-(\d+)", 1).alias("gene"),
                pl.col("name").str.extract(r"(.+)_(.+):(\d+)-(\d+)", 2).alias("transcript_ori"),
                pl.col("name").str.extract(r"(.+)_(.+):(\d+)-(\d+)", 3).cast(pl.UInt32).alias("pos_start"),
                pl.col("name").str.extract(r"(.+)_(.+):(\d+)-(\d+)", 4).cast(pl.UInt32).alias("pos_end"),
            ]
            if split_name
            else []
        )
        .with_columns(
            [
                (pl.col("transcript") == pl.col("transcript_ori")).alias("is_ori_seq"),
            ]
            + [
                (pl.col("pos_end") - pl.col("pos_start") + 1).alias("length"),
            ]
            if split_name
            else []
        )
        .collect()
    )
    return df


# %%
def stripplot(**kwargs: Any) -> Axes:
    import pandas as pd
    import seaborn as sns

    sns.set()

    df = pd.concat(pd.DataFrame({"x": v, "y": k}) for k, v in kwargs.items())
    return sns.stripplot(data=df, x="x", y="y", **(dict(orient="h", alpha=0.6)))
