# tblastn-lite

Find the DNA that encodes a protein. One C++ file, no dependencies, no database
build step, no BLAST install: point it at a protein FASTA and a nucleotide FASTA
and it does a 6-frame translated search.

820-aa protein vs the whole *E. coli* K-12 genome (4.6 Mb): **0.18 s** on 8 threads.

## Build

```bash
make
```

Or directly:

```bash
g++ -O3 -std=c++17 -static -o tblastn_lite tblastn_lite.cpp -pthread
```

`-static` gives a single portable `.exe` on Windows/MinGW.

## Use

```bash
./tblastn_lite -q protein.faa -d genome.fna -e 1e-5
```

```
# tblastn-lite  db: K12.fasta (1 seqs, 4641652 bp)
# query: thrA_ecoli (820 aa)  hits: 12
# qid        sid          %id   len  mism  frame  qstart qend sstart  send    evalue    bits
thrA_ecoli   NC_000913.3  100.0 820  0     +1     1      820  337     2796    0         1902.7
thrA_ecoli   NC_000913.3  47.5  158  83    +2     550    707  4131464 4131937 7.4e-46   182.7
thrA_ecoli   NC_000913.3  45.5  112  61    -3     184    295  4232702 4232367 3.8e-29   127.2
```

Subject coordinates are 1-based on the plus strand; `sstart > send` means the
match is on the minus strand.

Add `--aln` for the alignment and `--dna` to print the matching DNA as FASTA —
that second one is usually the point:

```bash
./tblastn_lite -q protein.faa -d genome.fna -e 1e-20 --dna > gene_candidates.fna
```

| flag | default | meaning |
|------|---------|---------|
| `-q, --query` | – | protein FASTA (multiple queries fine) |
| `-d, --db` | – | nucleotide FASTA (multiple sequences fine) |
| `-e, --evalue` | 10 | E-value cutoff |
| `-T, --thresh` | 13 | neighbourhood word threshold; **lower = more sensitive, slower** |
| `-X, --xdrop` | 16 | ungapped extension X-drop |
| `-t, --threads` | all cores | threads |
| `-m, --max` | 500 | max HSPs reported per query |
| `--aln` | off | print alignments |
| `--dna` | off | print matching DNA as FASTA |
| `--selftest` | – | built-in checks |

## How it works

Same recipe as the original BLAST, minus the parts you rarely need:

1. Index every 3-mer that scores ≥ `-T` against the query under BLOSUM62
   (the neighbourhood, not just exact words — that is where the sensitivity
   comes from).
2. Translate each subject sequence in all 6 frames, look every 3-mer up.
3. Extend each seed left and right without gaps until the score drops `-X`
   below its best; one HSP per diagonal.
4. Score with Karlin–Altschul (λ=0.318, K=0.134): `E = 2·m·n·2^-bits`.

Finds a 35%-substituted version of a query at E ≈ 1e-303, and its distant
paralogues at 1e-12.

## What it does not do

- **No gapped alignment.** HSPs are ungapped, so an indel splits a gene into
  two HSPs on adjacent diagonals. Fine for bacterial genomes and exact/near
  homologues; if you need gapped, composition-adjusted, frameshift-tolerant
  alignment, use NCBI `tblastn`.
- **No low-complexity (SEG) filter.** A query full of `PPPPQQQQ` will hit noise.
- **No spliced alignment.** Eukaryotic genes with introns give one HSP per exon.
- Sequences are held in memory (genome-sized is fine; don't feed it all of nt).

## Test

```bash
make test
```

`--selftest` checks translation, reverse complement, and an end-to-end search
for a protein planted in random DNA on both strands, asserting full-length
recovery and exact coordinate round-tripping.

## Licence

MIT.
