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


def nonspecific_test(probe: str, seq: str, t: float = 37):
    model = gen_model(t)
    probe_ = Strand(reverse_complement("TTT" + probe + "TTT"), "probe")
    seq_ = Strand(seq, "seq")
    t1 = Tube(
        strands={seq_: 1e-10, probe_: 1e-9}, name="t1", complexes=SetSpec(max_size=2)
    )
    # return tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    result: Result = tube_analysis(tubes=[t1], compute=["mfe"], model=model)
    tube_result: TubeResult = result[t1]
    want = {}
    for complex, conc in tube_result.complex_concentrations.items():
        if complex.name == "(probe)":
            want["hairpin"] = conc / 1e-9
        elif complex.name == "(probe+seq)" or complex.name == "(seq+probe)":
            want["bound"] = conc / 1e-10
        elif complex.name == "(probe+probe)":
            want["dimer"] = conc / 1e-9
    return result, want
