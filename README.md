# blastannotate

`blastannotate` is a high-performance C++ tool for performing BLASTn alignments, extracting top hits & annotating results with species information via NCBI Entrez Direct utilities. Designed for microbial genomics & comparative genomics pipelines, it facilitates reproducible & automated species annotation workflows.

## Features
- Runs `blastn` using a local nucleotide database
- Extracts top 10 hits per query
- Annotates hits with species name using NCBI Entrez utilities (`efetch`, `esearch`, `elink`, `xtract`)
- Generates human-readable summary reports (`TSV + summary`)
- **Multi-threaded** BLAST execution.  


## Requirements
- C++17 Compiler (`g++`/`clang++`) 
- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [NCBI Entrez Direct (EDirect)](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

## Installation

# blastannotate  

**Annotate BLAST results with species information & generate reports.**  

---

### 1. Clone & compile  

```bash

```bash
git clone https://github.com/gmboowa/blastannotate.git

## C++ installation
cd blastannotate/bin
make
g++ -std=c++17 -O3 -lpthread -o blastannotate blastannotate.cpp
sudo mv blastannotate /usr/local/bin/  # Optional: Install globally

```
### 2. Verify installation

```bash

blastannotate --help

```
## Usage 

```bash

blastannotate --ref_db your_db.fasta --output_dir output -fasta_1 query1.fa -fasta_2 query2.fa

Or bash usage

chmod 777 blastannotate

./blastannotate.sh --ref_db your_db.fasta --output_dir output -fasta_1 query1.fa -fasta_2 query2.fa

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



