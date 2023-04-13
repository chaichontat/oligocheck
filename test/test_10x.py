from oligocheck.tenx import gen_oligo
from oligocheck.oligocheck.sequtils import gen_random_base, pcr, reverse_complement

truncated_truseq_1 = "CTACACGACGCTCTTCCGATCT"
cell_bc = lambda: gen_random_base(16)
umi = lambda: gen_random_base(12)
polyA = lambda: "A" * 30
bead = lambda: truncated_truseq_1 + cell_bc() + umi()
rev_tso = reverse_complement("AAGCAGTGGTATCAACGCAGAGTACATGGG")


def gen_mrna():
    # https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3fb.html
    # http://tucf-genomics.tufts.edu/documents/protocols/TUCF_Understanding_Illumina_TruSeq_Adapters.pdf
    return bead() + polyA() + gen_random_base(200) + rev_tso


def test_mrna():
    mrna = gen_mrna()
    fwd = "CTACACGACGCTCTTCCGATCT"
    rev = "AAGCAGTGGTATCAACGCAGAG"
    assert pcr(mrna, fwd, rev)


def gen_adt():
    seq = gen_oligo([1], "polyA", "ADT").split("\t")[1][8:].replace("*", "")
    return bead() + reverse_complement(seq) + rev_tso


def test_adt():
    adt = gen_adt()
    fwd = "CTACACGACGCTCTTCCGATCT"
    rev = "CCTTGGCACCCGAGAATTCC"
    assert pcr(adt, fwd, rev)
