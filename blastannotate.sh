#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 --ref_db <DATABASE_FASTA> --output_dir <OUTPUT_DIR> -fasta_1 <QUERY_FASTA_1> [-fasta_2 <QUERY_FASTA_2> ... -fasta_n <QUERY_FASTA_N>]"
    echo "Example: $0 --ref_db /path/to/database.fasta --output_dir /path/to/output -fasta_1 /path/to/sample1.fa -fasta_2 /path/to/sample2.fa"
    exit 1
}

# Parse command-line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ref_db)
            DATABASE_FASTA="$2"
            shift 2
            ;;
        --output_dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -fasta_*)
            QUERY_FILES+=("$2")
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check for required arguments
if [ -z "$DATABASE_FASTA" ] || [ -z "$OUTPUT_DIR" ] || [ ${#QUERY_FILES[@]} -eq 0 ]; then
    echo "Error: Missing required arguments"
    usage
fi

DB_NAME=$(basename "$DATABASE_FASTA" | cut -d'.' -f1)

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Function to check BLAST indexes
check_blast_index() {
    local db="$1"
    local db_dir=$(dirname "$db")
    local db_base=$(basename "$db")
    
    local nhr_count=$(find "$db_dir" -name "${db_base}*.nhr" | wc -l)
    local nin_count=$(find "$db_dir" -name "${db_base}*.nin" | wc -l)
    local nsq_count=$(find "$db_dir" -name "${db_base}*.nsq" | wc -l)
    
    if [[ $nhr_count -gt 0 && $nin_count -gt 0 && $nsq_count -gt 0 ]]; then
        blastdbcmd -db "$db" -info >/dev/null 2>&1 && return 0
    fi
    return 1
}

# Step 1: Setup BLAST database
if check_blast_index "$DATABASE_FASTA"; then
    echo "BLAST database is properly indexed - skipping indexing"
else
    echo "Indexing BLAST database..."
    makeblastdb -in "$DATABASE_FASTA" -dbtype nucl -parse_seqids \
               -title "${DB_NAME}_DB" -max_file_sz "2GB" || exit 1
fi

# Function to extract clean accession
clean_accession() {
    echo "$1" | sed -e 's/^[^|]*|//' -e 's/|.*//' -e 's/\..*//'
}

# Function to get species from accession
get_species() {
    local acc="$1"
    local species=""
    local max_retries=3
    local retry_delay=1
    local clean_acc=$(clean_accession "$acc")
    
    for ((i=1; i<=$max_retries; i++)); do
        # First try direct nucleotide record
        species=$(efetch -db nuccore -id "$clean_acc" -format docsum 2>/dev/null | \
                 xtract -pattern DocumentSummary -element Organism 2>/dev/null | \
                 head -1)
        
        # If that fails, try taxonomy database
        if [ -z "$species" ]; then
            taxid=$(esearch -db nuccore -query "$clean_acc" 2>/dev/null | \
                   elink -target taxonomy 2>/dev/null | \
                   efetch -format uid 2>/dev/null | \
                   head -1)
            
            if [ -n "$taxid" ]; then
                species=$(efetch -db taxonomy -id "$taxid" -format docsum 2>/dev/null | \
                         xtract -pattern DocumentSummary -element ScientificName 2>/dev/null)
            fi
        fi
        
        [ -n "$species" ] && break
        sleep $retry_delay
    done
    
    [ -z "$species" ] && species="Unknown"
    echo "$species"
}

# Process each query file
for QUERY_FASTA in "${QUERY_FILES[@]}"; do
    QUERY_NAME=$(basename "$QUERY_FASTA" | cut -d'.' -f1)
    SAMPLE_DIR="${OUTPUT_DIR}/${QUERY_NAME}"
    mkdir -p "$SAMPLE_DIR"
    
    echo -e "\nProcessing sample: $QUERY_NAME"
    
    # Run BLASTn
    BLAST_OUTPUT="${SAMPLE_DIR}/blast_results.tsv"
    echo " - Running BLASTn..."
    blastn -query "$QUERY_FASTA" -db "$DATABASE_FASTA" \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
           -num_threads $(nproc) -out "$BLAST_OUTPUT" || continue
    
    # Get top 10 hits
    TOP_HITS="${SAMPLE_DIR}/top10_hits.tsv"
    echo " - Extracting top 10 hits..."
    (echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore"; 
     sort -k12,12gr "$BLAST_OUTPUT" | head -10) > "$TOP_HITS"
    
    # Annotate with species
    ANNOTATED_HITS="${SAMPLE_DIR}/top10_hits_with_species.tsv"
    echo " - Annotating with species..."
    
    # Create species mapping
    cut -f2 "$TOP_HITS" | tail -n +2 | while read -r acc; do
        clean_acc=$(clean_accession "$acc")
        species=$(get_species "$clean_acc")
        echo -e "$acc\t$species"
        sleep 0.5
    done > "$SAMPLE_DIR/species_mapping.tsv"
    
    # Merge with top hits
    awk -F'\t' '
    BEGIN {OFS="\t"}
    NR==FNR {map[$1]=$2; next}
    FNR==1 {print $0, "species_name"; next}
    {
        print $0, (map[$2] ? map[$2] : "Unknown")
    }' "$SAMPLE_DIR/species_mapping.tsv" "$TOP_HITS" > "$ANNOTATED_HITS"
    
    # Generate report
    REPORT_FILE="${SAMPLE_DIR}/report.txt"
    echo " - Generating report..."
    (
        echo "BLAST Analysis Report"
        echo "===================="
        echo "Sample: $QUERY_NAME"
        echo "Date: $(date)"
        echo "Database: $DB_NAME"
        echo ""
        echo "Top Hits with Species:"
        column -t -s $'\t' "$ANNOTATED_HITS"
    ) > "$REPORT_FILE"
    
    echo " - Results saved to $SAMPLE_DIR"
done

echo -e "\nAll samples processed successfully!"
echo "Final output directories:"
find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
    echo "- $dir"
done

exit 0