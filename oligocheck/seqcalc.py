import io
import re
from functools import cache

import polars as pl
import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

from oligocheck.merfish.pairwise import pairwise_alignment
from oligocheck.sequtils import gc_content

deltas = pl.read_csv(
    io.StringIO(
        """
key,G,H,S
rAAIdTT,-0.89,-6.62,-18.5
rAUIdTA,-0.89,-6.37,-17.7
rAGIdTC,-1.52,-8.30,-21.8
rACIdTG,-1.77,-9.79,-25.8
rUA/dAT,-0.60,-9.89,-30.0
rUU/dAA,-0.34,-5.33,-16.1
rUG/dAC,-1.39,-5.33,-12.7
rUC/dAG,-1.17,-9.29,-26.1
rGA/dCT,-1.69,-10.78,-29.3
rGU/dCA,-1.51,-9.26,-25.0
rGG/dCC,-2.52,-9.03,-20.9
rGC/dCG,-2.54,-14.06,-37.2
rCA/dGT,-1.36,-7.15,-18.7
rCU/dGA,-0.99,-4.79,-12.2
rCG/dGC,-1.67,-8.53,-22.1
rCC/dGG,-1.87,-5.28,-11.0
init,+2.65,-3.11,-18.8
""".replace(
            "I", "/"
        )
        .replace("U", "T")
        .replace("r", "")
        .replace("d", "")
    )
)

R_DNA = dict(zip(deltas["key"], deltas[["H", "S"]].iter_rows())) | {
    "init_A/T": (0, 0),
    "init_G/C": (0, 0),
    "init_oneG/C": (0, 0),
    "init_allA/T": (0, 0),
    "init_5T/A": (0, 0),
    "sym": (0, 0),
}

CONDITIONS = {
    "fish": dict(nn_table=mt.DNA_NN3, Na=330, Tris=0, Mg=0, dNTPs=0, dnac1=1, dnac2=0),
    "hybrid": dict(nn_table=R_DNA, Na=330, Tris=0, Mg=0, dNTPs=0, dnac1=1, dnac2=0),
    "q5": dict(nn_table=mt.DNA_NN3, Na=100, Tris=0, Mg=2, dNTPs=0.2, dnac1=500, dnac2=10),
}


def tm_q5(seq: str, **kwargs: float) -> float:
    return mt.Tm_NN(Seq(seq), **{**CONDITIONS["q5"], **kwargs})


def formamide_molar(percent: float = 30) -> float:
    return percent * 10 * 1.13 / 45.04  # density / molar mass


def formamide_correction(seq: str, fmd: float = 30) -> float:
    return (0.453 * (gc_content(seq) / 100.0) - 2.88) * formamide_molar(fmd)


def tm_fish(seq: str, **kwargs: float) -> float:
    return mt.Tm_NN(Seq(seq), **{**CONDITIONS["fish"], **kwargs}) + formamide_correction(seq)


def tm_hybrid(seq: str, **kwargs: float) -> float:
    return mt.Tm_NN(Seq(seq), **{**CONDITIONS["hybrid"], **kwargs}) + formamide_correction(seq)


def hp_fish(seq: str, **kwargs: float) -> float:
    return primer3.calc_hairpin_tm(seq, **{**dict(mv_conc=330, dv_conc=0, dntp_conc=0, dna_conc=1), **kwargs})


r = re.compile(r"(\|{15,})")


@cache
def tm_match(seq: str, cigar: str, mismatched_reference: str) -> float:
    p = pairwise_alignment(seq, cigar, mismatched_reference)
    try:
        return max(tm_hybrid(seq[slice(*x.span())]) for x in r.finditer(p[1]))
    except ValueError:
        return 0
