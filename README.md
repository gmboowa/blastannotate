# blastannotate

`blastannotate` is a high-performance C++ tool for performing BLASTn alignments, extracting top hits & annotating results with species information via NCBI Entrez Direct utilities. Designed for microbial genomics & comparative genomics pipelines, it facilitates reproducible & automated species annotation workflows.

## Features
- Runs `blastn` using a local nucleotide database
- Extracts top 10 hits per query
- Annotates hits with species name using NCBI Entrez utilities (`efetch`, `esearch`, `elink`, `xtract`)
- Generates human-readable summary reports

## Requirements
- `g++` (C++17 or later)
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [NCBI Entrez Direct (EDirect)](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

## Installation

```bash
git clone https://github.com/gmboowa/blastannotate.git
cd blastannotate
make
```

## Usage

```bash
./blastannotate --ref_db your_db.fasta --output_dir output \
    -fasta_1 query1.fa -fasta_2 query2.fa
```

## Output

For each sample:

```
output/
├── sample1/
│   ├── blast_results.tsv
│   ├── top10_hits_with_species.tsv
│   └── report.txt
└── sample2/
    └── ...
```



