from typing import Callable, Literal

from .io import get_whitelist

whitelist = get_whitelist()
# https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3fb.html
# https://www.biolegend.com/en-us/products/totalseq-b0162-anti-human-cd64-antibody-18415
# https://www.biolegend.com/en-us/products/totalseq-b0302-anti-mouse-hashtag-2-antibody-17772
FLANKING_5 = {
    # Same 5' labeling as all other labeling, differentiation is in the 3' capture sequence.
    "TruSeq Read 2": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNN",
    # Small RNA - ADT
    "RPI1": "CCTTGGCACCCGAGAATTCCA",
    # TruSeq D701_s - Hashtag
    "D701": "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
    "NexteraTranposase": "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
}

FLANKING_3 = {
    "polyA": "BAAAAAAAAAAAAAAAAAAAAAAAAAAAA*A*A",
    "Capture1": "NNNNNNNNNGCTTTAAGGCCGGTCCTAGC*A*A",
    "Capture2": "GCTCACCTATTAGCGGCTAA*G*G",
}


def gen_flank(
    provider: Literal["10x", "polyA", "sci"], type_: Literal["ADT", "HTO"] | None = None
) -> tuple[str, str]:
    if provider == "10x":
        if type_ == "ADT":
            return FLANKING_5["TruSeq Read 2"], FLANKING_3["Capture1"]
        elif type_ == "HTO":
            return FLANKING_5["TruSeq Read 2"], FLANKING_3["Capture2"]
        else:
            raise ValueError(f"Unknown type: {type_}")

    elif provider == "polyA":
        if type_ == "ADT":
            return FLANKING_5["RPI1"], FLANKING_3["polyA"]
        elif type_ == "HTO":
            return FLANKING_5["D701"], FLANKING_3["polyA"]
        else:
            raise ValueError(f"Unknown type: {type_}")
    elif provider == "sci":
        return FLANKING_5["NexteraTranposase"], FLANKING_3["polyA"]
    else:
        raise ValueError(f"Unknown provider: {provider}")


# https://www.biolegend.com/en-us/protocols/totalseq-a-antibodies-and-cell-hashing-with-10x-single-cell-3-reagent-kit-v3-3-1-protocol
# TotalSeq-A
# -------Flank--------|------BC------|------------polyA-----------------|
# CCTTGGCACCCGAGAATTCCAAACAAGACCCTTGAGBAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*A*A

# Illumina Small RNA RPI1 primer
# (for ADT amplification; i7 index 1, Oligonucleotide sequences, Illumina)
# -----------------------|Index|----------------------------------|
# CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTC*C*A

# Illumina TruSeq D701_LONG primer
# (for HTO amplification; i7 index 1)
# CAAGCAGAAGACGGCATACGAGATCGAGTAATGTGACTGGAGTTCAGACGTGTGCTCTTCCGAT*C*T

# sci-rna-seq3
# Nextera Transposase Adapters Read 2
# GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG[10bp-barcode]BAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA*A*A


def default_builder(
    seq: str,
    i: int,
    modification: str = "/5AmMC12/",
) -> str:
    return "\t".join([modification + seq])


def default_namer(x: str, i: int) -> str:
    return f"{i}"


def gen_oligo(
    idxs: list[int],
    *,
    builder: Callable[[str, int], str] = default_builder,
    namer: Callable[[str, int], str] = default_namer,
    provider: Literal["10x", "polyA", "sci"] = "10x",
    scale: str = "100nm",
    purification: str = "STD",
    type_: Literal["ADT", "HTO"] | None = None,
) -> str:
    f5, f3 = gen_flank(provider, type_)
    if 0 in idxs:
        raise ValueError("This is 1-indexed.")

    seqs = ["\t".join([builder(f"{f5}{whitelist[i-1]}{f3}", i), scale, purification]) for i in idxs]
    names = [namer(whitelist[i - 1], i) for i in idxs]
    return "\n".join(["\t".join(x) for x in zip(names, seqs)])
