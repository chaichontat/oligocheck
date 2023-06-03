import logging
from functools import cache

from oligocheck.sequtils import reverse_complement

logger = logging.getLogger("nupack")
logger.propagate = False

from nupack import Model, SetSpec, Strand, Tube, tube_analysis  # noqa: E402
from nupack.analysis import Result, TubeResult  # noqa: E402


@cache
def gen_model(
    t: float,
    formamide: float = 30,
    sodium: float = 0.3,
    magnesium: float = 0.0,
    **kwargs,
):
    return Model(
        material="dna",
        celsius=t + formamide * 0.65,
        sodium=sodium,
        magnesium=magnesium,
        **kwargs,
    )


def nonspecific_test(
    probe: str,
    seq: str,
    t: float = 37,
    formamide: float = 30,
    conc_seq: float = 1e-10,
    conc_probe: float = 1e-9,
):
    model = gen_model(t, formamide=formamide)
    probe_ = Strand(reverse_complement(probe), "probe")
    seq_ = Strand(seq, "seq")
    t1 = Tube(strands={seq_: conc_seq, probe_: conc_probe}, name="t1", complexes=SetSpec(max_size=2))
    # return tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    result: Result = tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    tube_result: TubeResult = result[t1]
    want = {}
    for complex, conc in tube_result.complex_concentrations.items():
        if complex.name == "(probe)":
            want["hairpin"] = conc / conc_probe
        elif complex.name == "(probe+seq)" or complex.name == "(seq+probe)":
            want["bound"] = conc / conc_seq
        elif complex.name == "(probe+probe)":
            want["dimer"] = conc / conc_probe
    return result, want


def secondary_structure(
    probe: str,
    t: float = 37,
    formamide: float = 30,
    conc_probe: float = 1e-9,
):
    model = gen_model(t, formamide=formamide)
    probe_ = Strand(reverse_complement(probe), "probe")
    t1 = Tube(strands={probe_: conc_probe}, name="t1", complexes=SetSpec(max_size=2))
    # return tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    result: Result = tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    return result
