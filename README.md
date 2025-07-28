# blastannotate

`blastannotate` is a bash tool for performing BLASTn alignments, extracting top hits & annotating results with species information via NCBI Entrez Direct utilities. Designed for microbial genomics & comparative genomics pipelines, it facilitates reproducible & automated species annotation workflows.

## Features
- Runs `blastn` using a local nucleotide database
- Extracts top 10 hits per query
- Annotates hits with species name using NCBI Entrez utilities (`efetch`, `esearch`, `elink`, `xtract`)
- Generates human-readable summary reports (`TSV + summary`)
- **Multi-threaded** BLAST execution.  


## Requirements

- [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
- [NCBI Entrez Direct (EDirect)](https://www.ncbi.nlm.nih.gov/books/NBK179288/)

## Installation

**Annotate BLAST results with species information & generate reports.**  

---

### 1. Clone & compile  

```bash

# Clone the repository
git clone https://github.com/gmboowa/blastannotate.git

# Navigate into the project directory
cd blastannotate

# Make the shell script executable (if using a wrapper script)
chmod +x blastannotate.sh


```
### 2. Verify installation

```bash


./blastannotate.sh --help

```
## Usage 

```bash


./blastannotate.sh --ref_db your_db.fasta --output_dir output -fasta_1 query1.fa -fasta_2 query2.fa

```

## Output

Sample of the outputs:

```bash
output/
NCBI tools missing: efetch esearch elink xtract
Attempting to install...
Installing NCBI E-utilities...
Found existing edirect installation, adding to PATH...
BLAST database is properly indexed - skipping indexing

Processing sample: A55727
 - Running BLASTn...
 - Extracting top 10 hits...
 - Annotating with species...
 - Generating report...
 - Results saved to /A55727

Processing sample: A56429
 - Running BLASTn...
 - Extracting top 10 hits...
 - Annotating with species...
 - Generating report...
 - Results saved to /A56429

Processing sample: ERR2870192
 - Running BLASTn...
 - Extracting top 10 hits...
 - Annotating with species...
 - Generating report...
 - Results saved to /ERR2870192

Processing sample: ERR13420967
 - Running BLASTn...
 - Extracting top 10 hits...
 - Annotating with species...
 - Generating report...
 - Results saved to /ERR13420967

Processing complete!
Output directories created:
- /ERR2870192
- A56429
- /A55727
- ERR13420967

```

## Sample BLAST Report Output

```bash
BLAST Analysis Report
====================
Sample: A55727
Date: 2025-07-28 12:17:14
Database: eskapee_combined

Top Hits with Species:

qseqid         sseqid              pident    length    mismatch  gapopen   qstart    qend      sstart    send      evalue         bitscore       species_name                  
k121_326       ref|NZ_CP014071.1|  99.998    116710    2         0         1         116710    264899    381608    0.00e+00       2.2e+05        Klebsiella quasipneumoniae    
k121_326       ref|NZ_CP029597.1|  1e+02     116710    2         0         1         116710    2720093   2603384   0.00e+00       2.2e+05        Klebsiella quasipneumoniae subsp. similipneumoniae
k121_326       ref|NZ_CP014696.2|  1e+02     116710    2         0         1         116710    3198456   3315165   0.00e+00       2.2e+05        Klebsiella quasipneumoniae    
k121_326       ref|NZ_CP084777.1|  1e+02     116712    124       5         1         116710    573345    456643    0.00e+00       2.1e+05        Klebsiella quasipneumoniae subsp. similipneumoniae
k121_326       ref|NZ_CP077279.1|  1e+02     116712    124       5         1         116710    4604222   4720924   0.00e+00       2.1e+05        Klebsiella quasipneumoniae    
k121_326       ref|NZ_CP047283.1|  1e+02     116712    125       5         1         116710    572971    456269    0.00e+00       2.1e+05        Klebsiella quasipneumoniae    
k121_326       ref|NZ_CP063896.1|  1e+02     116701    427       7         14        116710    3017055   2900364   0.00e+00       2.1e+05        Klebsiella quasipneumoniae    
k121_326       ref|NZ_CP152902.1|  1e+02     116701    419       5         14        116710    577898    461205    0.00e+00       2.1e+05        Klebsiella quasipneumoniae subsp. similipneumoniae
k121_326       ref|NZ_CP140611.1|  1e+02     116698    425       4         14        116710    567972    451284    0.00e+00       2.1e+05        Klebsiella quasipneumoniae    
k121_326       ref|NZ_CP031257.1|  1e+02     116698    427       4         14        116710    569817    453129    0.00e+00       2.1e+05        Klebsiella quasipneumoniae    
---

