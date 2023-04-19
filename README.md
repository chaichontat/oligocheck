# oligocheck

A barcode designed and PCR simulator

# MERFISH Probes

The reference is Ensembl release 109 mm39 [cDNA](https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/) and [ncRNA](https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/ncrna/).
tRNA information is from [GtRNAdb](http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc39/Mmusc39-seq.html) high confidence data release 20.
All FASTA files were concatenated and passed into bowtie2 for index generation.

```sh
bowtie2-build combi.fa mm39
```

## Generate codebook
The MHD4 code is generated using the MATLAB script from the authors.
The output is a CSV table.
The FPKM count from the mouse embryonic brain [ENCODE](https://www.encodeproject.org/publication-data/ENCSR574CRQ/) data is used to assign each gene to each barcode.
Put simply, shuffle, sum up the total expression of each bit, and minimize the entropy.


## Get probe list

The input is a list of gene symbols. As of now, we do not care about isoforms,
so I decided to use the [GENCODE basic](https://www.gencodegenes.org/mouse/) annotation release M32 to pick the most common and least suspect isoforms.
The GTF file is filtered to retain only transcript information.
The input list is queried through a Ensembl REST API to get the list of ensembl gene IDs, which is used to look up transcript IDs in the filtered database.
The sequences are obtained using the Ensembl REST API (although we should probably use `pyfaidx`) going forward.
All lookups are cached locally using Redis.

The sequences are passed into `blockParse.py` in overlap mode, which slides through the entire sequence and extracts every substring that is within a set Tm (DNA_NN3) 49-57°C with 30% formamide.
This is to maximize the number of potential candidates which can be filtered out later.
Probe length is between 25-41 bp.
These are passed into `bowtie2` with arguments

```python
# --no-hd No SAM header
# -k 100 report first 100 alignments
# -D 20 consecutive seed extension attempts can "fail" before Bowtie 2 moves on
# -R 3 the maximum number of times Bowtie 2 will "re-seed" reads with repetitive seeds.
# -L 17 seed length
# -i C,2 Seed interval, every 2 bp
# --score-min G,1,4 f(x) = 1 + 4*ln(read_length)
f"bowtie2 -x data/mm39/mm39 {(temp / f'{ts}.fastq').as_posix()} "
f"--no-hd -k 100 --local -D 20 -R 3 -N 0 -L 17 -i C,2 --score-min G,1,4 "
f"-S {(temp / f'{ts}.sam').as_posix()}",
```

Effectively, this is 'detect anything with homology greater than or equal to 17 bp'.
The output SAM file is parsed. Pool all probes from the same gene across isoforms and filter using the mapped transcript identity:

- If from the same gene, count if within GENCODE basic, otherwise ignore. Assuming that it'll bind to the isofrom.
- If mapped to rRNA or tRNA, discard right away.

Next, filter based on sequences, the first two is to save some time, since these often have unacceptably high Tms regardless.

- If the Levenshtein distance between the substring of the mapped transcript and the probe is less than 5, discard.
- If the max of CIGAR matches is greater than or equal to 20, discard.
- Perform NUPACK simulation of binding at 47°C + formamide correction. If steady-state concentration > 5% of total target, discard.

Then, score each probe with desirable characteristics from [https://www.pnas.org/doi/10.1073/pnas.0812506106]:

> On the basis of these nucleotide composition analyses, we derived 2 new probe design rules:
> (i) to improve probe responsiveness, the nucleotide composition of A in a probe should be limited to below 28%, and AAAA stacks should be avoided in probe sequences;
> (ii) to reduce cross-hybridization effects but still maintain reasonable probe response, the C nucleotide composition of probes should be limited to between 22 and 28%,
> and CCCC stack or 4 nonconsecutive Cs in any 6 consecutive nucleotides in the first 12 positions of a probe should be avoided.

and GC content between 35-65%.

## Construct encoding probes

All pairs of readout sequence are screened against the transcriptome using bowtie2 as described above.
Using CIGAR to extract matches, anything with Tm above 33 is discarded.

For each probe, cycle through the permutation of its readout sequence, and concatenate with readout sequences with TT spacing (to be reverse complemented). Then concatenate with PCR header and the necessary [spacers](https://www.nature.com/articles/s42003-020-01167-x) for the T7 promoter

Pass the entire thing into bowtie2 again.

Off-target binding with Tm >37°C is discarded.

Sort the entire probe sequence by scores from (off-target binding with Tm > 33°C), score from the previous section, and number of isoforms mapped, and pick the top 40 probes.
Concatenate with the remaining of the T7 promoter and PCR handle.
