import primer3
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt

CONDITIONS = {
    "fish": dict(nn_table=mt.DNA_NN3, Na=330, Tris=0, Mg=0, dNTPs=0, dnac1=1, dnac2=0),
    "q5": dict(nn_table=mt.DNA_NN3, Na=100, Tris=0, Mg=2, dNTPs=0.2, dnac1=500, dnac2=10),
}


def tm_q5(seq: str, **kwargs: float) -> float:
    return mt.Tm_NN(Seq(seq), **{**CONDITIONS["q5"], **kwargs})


def tm_fish(seq: str, **kwargs: float) -> float:
    return mt.Tm_NN(Seq(seq), **{**CONDITIONS["fish"], **kwargs})


def hp_fish(seq: str, **kwargs: float) -> float:
    return primer3.calc_hairpin_tm(seq, **{**dict(mv_conc=330, dv_conc=0, dntp_conc=0, dna_conc=1), **kwargs})
