#%%
import random
from binascii import hexlify
from collections import Counter
from itertools import product
from typing import Callable

import colorama
import pandas as pd
import primer3
import seaborn as sns
from hexhamming import hamming_distance_string
from Levenshtein import distance
from nupack import *

#%%
wells = [f"{row}{col}" for row in "ABCDEFGH" for col in range(1, 13)]


def gen_plate(name: str, seqs: list[str]):
    wells = [f"{row}{col}" for row in "ABCDEFGH" for col in range(1, 13)]
    return pd.DataFrame(
        {
            "Well Position": wells,
            "Name": [name + w for w in wells],
            "Sequence": seqs,
        }
    )


def reverse_complement(seq: str):
    return seq[::-1].translate(str.maketrans("ATCG", "TAGC"))


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


def gc_content(seq: str):
    return (seq.count("G") + seq.count("C")) / len(seq)


def calculate_base_diversity(seq: list[str]):
    return [{k: v / len(seq) for k, v in Counter([x[i] for x in seq]).items()} for i in range(len(seq[0]))]


def gen_sample_space(template: str, seed: int = 42):
    mapping = {
        "A": "A",
        "T": "T",
        "G": "G",
        "C": "C",
        "N": "ATGC",
        "W": "AT",
        "S": "GC",
        "M": "AC",
        "K": "GT",
        "R": "AG",
        "Y": "CT",
        "B": "CGT",
        "D": "AGT",
        "H": "ACT",
        "V": "ACG",
    }
    iters = [mapping[x] for x in template]
    base = ["".join(x) for x in product(*iters)]
    random.seed(42)
    random.shuffle(base)
    return base


def screen_umi_interactions(
    seqs: list[str],
    umis: list[str],
    conditions: dict[str, float],
    f: Callable[[str, str], str] = lambda x, y: x,
) -> list[list[float]]:
    idxs = []
    for seq in seqs:
        o = [primer3.calc_hairpin_tm(f(seq, umi), **conditions) for umi in umis]
        idxs.append(o)
    return idxs


def check_min_distance(seqs: list[str]):
    out = float("inf")
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            out = min(out, distance(seqs[i], seqs[j]))
    return out


def gen_rt():
    bridges = pd.read_csv("./ps_bridges.tsv", sep="\t")
    htvn = "".join(["H"] * 12) + "".join(["T"] * 23) + "VN"
    rtplate = gen_plate(
        "sciv2-RT-", ["/5Phos/CTCACTG" + bridges.iloc[s[1], 1][:10].lower() + htvn for s in selected]
    )
