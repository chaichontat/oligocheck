from oligocheck.utils import gen_random_base, pcr, reverse_complement

cell_bc = lambda: gen_random_base(16)
umi = lambda: gen_random_base(12)
polyA = lambda: "A" * 30


def gen_mrna():
    # https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3fb.html
    # http://tucf-genomics.tufts.edu/documents/protocols/TUCF_Understanding_Illumina_TruSeq_Adapters.pdf
    truncated_truseq_1 = "CTACACGACGCTCTTCCGATCT"
    rev_tso = reverse_complement("AAGCAGTGGTATCAACGCAGAGTACATGGG")
    return truncated_truseq_1 + cell_bc() + umi() + polyA() + gen_random_base(200) + rev_tso


def test_mrna():
    mrna = gen_mrna()
    fwd = "CTACACGACGCTCTTCCGATCT"
    rev = "AAGCAGTGGTATCAACGCAGAG"
    assert pcr(mrna, fwd, rev)
