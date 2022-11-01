import pytest
from oligocheck.utils import pcr, reverse_complement


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"


def test_pcr():
    seq = "GCAGCGTCAGATGTGTATAAGAGACAGTCAGTTTTTGCTAGGACCGGCCTTAAAGC"
    primer = "GCAGCGTCAGATGTGTATAAGAGACAG"
    primer_rc = "GCTTTAAGGCCGGTCCTAGC"
    assert pcr(seq, primer, primer_rc) == seq
    with pytest.raises(ValueError):
        pcr(seq, primer, primer_rc + "T")

    assert pcr(seq, primer, primer_rc[1:]) == seq[:-1]
