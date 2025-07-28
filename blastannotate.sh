#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 --ref_db <DATABASE_FASTA> --output_dir <OUTPUT_DIR> -fasta_1 <QUERY_FASTA_1> [-fasta_2 <QUERY_FASTA_2> ... -fasta_n <QUERY_FASTA_N>]"
    echo "Example: $0 --ref_db /path/to/database.fasta --output_dir /path/to/output -fasta_1 /path/to/sample1.fa -fasta_2 /path/to/sample2.fa"
    exit 1
}

# Function to check and add to PATH if not present
add_to_path() {
    local dir="$1"
    if [[ ":$PATH:" != *":$dir:"* ]]; then
        export PATH="$dir:$PATH"
        echo "export PATH=\"$dir:\$PATH\"" >> ~/.zshrc
        echo "export PATH=\"$dir:\$PATH\"" >> ~/.bash_profile
    fi
}

# Function to install Homebrew if missing
install_homebrew() {
    if ! command -v brew &> /dev/null; then
        echo "Installing Homebrew..."
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.zshrc
        echo 'eval "$(/opt/homebrew/bin/brew shellenv)"' >> ~/.bash_profile
        eval "$(/opt/homebrew/bin/brew shellenv)"
    fi
}

# Function to install BLAST+
install_blast() {
    if ! command -v blastn &> /dev/null; then
        echo "Installing BLAST+..."
        install_homebrew
        brew install blast
        add_to_path "/opt/homebrew/bin"
    fi
}

# Function to install NCBI E-utilities
install_edirect() {
    local edirect_path="$HOME/edirect"
   
    if ! command -v efetch &> /dev/null; then
        echo "Installing NCBI E-utilities..."
       
        if [ -d "$edirect_path" ]; then
            echo "Found existing edirect installation, adding to PATH..."
            add_to_path "$edirect_path"
            return
        fi
       
        sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
       
        if [ ! -d "$edirect_path" ]; then
            echo "Error: edirect installation failed!"
            return 1
        fi
       
        echo "Configuring PATH for edirect..."
        echo "source ~/.bash_profile" >> ~/.bashrc
        echo "export PATH=\"$edirect_path:\$PATH\"" >> ~/.bash_profile
        echo "export PATH=\"$edirect_path:\$PATH\"" >> ~/.zshrc
        add_to_path "$edirect_path"
       
        for tool in efetch esearch elink xtract; do
            if ! command -v $tool &> /dev/null; then
                echo "Warning: $tool still not found after installation"
            fi
        done
    fi
}

# Function to verify NCBI tools are available
check_ncbi_tools() {
    local missing_tools=()
    for tool in efetch esearch elink xtract; do
        if ! command -v $tool &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
   
    if [ ${#missing_tools[@]} -gt 0 ]; then
        echo "NCBI tools missing: ${missing_tools[*]}"
        echo "Attempting to install..."
        install_edirect
       
        missing_tools=()
        for tool in efetch esearch elink xtract; do
            if ! command -v $tool &> /dev/null; then
                missing_tools+=("$tool")
            fi
        done
       
        if [ ${#missing_tools[@]} -gt 0 ]; then
            echo "Warning: Missing NCBI tools after installation: ${missing_tools[*]}"
            echo "Species annotation will be skipped"
            return 1
        fi
    fi
    return 0
}

# Function to verify all tools are available
verify_dependencies() {
    install_homebrew
    install_blast
    check_ncbi_tools
   
    local missing_tools=()
    for tool in blastn makeblastdb; do
        if ! command -v $tool &> /dev/null; then
            missing_tools+=("$tool")
        fi
    done
   
    if [ ${#missing_tools[@]} -gt 0 ]; then
        echo "Error: Critical tools missing: ${missing_tools[*]}"
        exit 1
    fi
}

# Function to clean accession numbers
clean_accession() {
    echo "$1" | sed -E 's/^[^|]*\|([^|]+)\|.*/\1/' | sed 's/\..*//'
}

# Function to get species from accession
get_species() {
    local acc="$1"
    local max_retries=3
    local retry_delay=2
    local clean_acc=$(clean_accession "$acc")
    local species="Unknown"
   
    for ((i=1; i<=$max_retries; i++)); do
        species=$(efetch -db nuccore -id "$clean_acc" -format docsum 2>/dev/null | \
                 xtract -pattern DocumentSummary -element Organism 2>/dev/null | \
                 head -1 | tr -d '\n')
       
        [ -n "$species" ] && break
       
        sleep $retry_delay
    done
   
    [ -z "$species" ] && species="Unknown"
    echo "$species"
}

# Function to batch process accessions (compatible version)
batch_get_species() {
    local input_file="$1"
    local output_file="$2"
    local batch_size=20
    local delay=3
    local count=0
   
    > "$output_file"  # Clear output file first
   
    while IFS= read -r acc; do
        {
            clean_acc=$(clean_accession "$acc")
            species=$(get_species "$clean_acc")
            echo -e "$acc\t$species"
        } >> "$output_file"
       
        ((count++))
        if (( count % batch_size == 0 )); then
            sleep $delay
        fi
    done < "$input_file"
}

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

# Enhanced argument parsing
declare -a QUERY_FILES
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
        *.fa|*.fasta|*.fna)
            if [[ "${#QUERY_FILES[@]}" -eq 0 ]]; then
                QUERY_FILES+=("$1")
                echo "Warning: Assuming $1 as -fasta_1 (recommend explicit -fasta_N flags)"
                shift
            else
                next_index=$(( ${#QUERY_FILES[@]} + 1 ))
                QUERY_FILES+=("$1")
                echo "Warning: Assuming $1 as -fasta_${next_index} (recommend explicit -fasta_N flags)"
                shift
            fi
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate arguments
if [ -z "$DATABASE_FASTA" ] || [ -z "$OUTPUT_DIR" ] || [ ${#QUERY_FILES[@]} -eq 0 ]; then
    echo "Error: Missing required arguments"
    usage
fi

# Verify and install all dependencies
verify_dependencies

# Check if NCBI tools are available for annotation
NCBI_TOOLS_AVAILABLE=true
if ! command -v efetch &> /dev/null || ! command -v xtract &> /dev/null; then
    NCBI_TOOLS_AVAILABLE=false
    echo "Warning: NCBI annotation tools not available - skipping species annotation"
fi

DB_NAME=$(basename "$DATABASE_FASTA" | cut -d'.' -f1)
mkdir -p "$OUTPUT_DIR"

# Setup BLAST database
if check_blast_index "$DATABASE_FASTA"; then
    echo "BLAST database is properly indexed - skipping indexing"
else
    echo "Indexing BLAST database..."
    makeblastdb -in "$DATABASE_FASTA" -dbtype nucl -parse_seqids \
               -title "${DB_NAME}_DB" -max_file_sz "2GB" || {
        echo "Error: Failed to index database"
        exit 1
    }
fi

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
           -num_threads $(nproc) -out "$BLAST_OUTPUT" || {
        echo "Warning: BLASTn failed for $QUERY_FASTA"
        continue
    }
   
    # Get top 10 hits
    TOP_HITS="${SAMPLE_DIR}/top10_hits.tsv"
    echo " - Extracting top 10 hits..."
    (echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore";
     sort -k12,12gr "$BLAST_OUTPUT" | head -10) > "$TOP_HITS"
   
    # Annotate with species if NCBI tools are available
    ANNOTATED_HITS="${SAMPLE_DIR}/top10_hits_with_species.tsv"
    if $NCBI_TOOLS_AVAILABLE; then
        echo " - Annotating with species..."
       
        # Create species mapping
        cut -f2 "$TOP_HITS" | tail -n +2 | sort -u > "$SAMPLE_DIR/unique_accessions.txt"
        batch_get_species "$SAMPLE_DIR/unique_accessions.txt" "$SAMPLE_DIR/species_mapping.tsv"
       
        # Merge with top hits
        if [ -f "$SAMPLE_DIR/species_mapping.tsv" ] && [ -s "$SAMPLE_DIR/species_mapping.tsv" ]; then
            awk -F'\t' '
            BEGIN {OFS="\t"}
            NR==FNR {map[$1]=$2; next}
            FNR==1 {print $0, "species"; next}
            {
                print $0, (map[$2] ? map[$2] : "Unknown")
            }' "$SAMPLE_DIR/species_mapping.tsv" "$TOP_HITS" > "$ANNOTATED_HITS"
        else
            cp "$TOP_HITS" "$ANNOTATED_HITS"
            sed -i.bak '1s/$/ species/' "$ANNOTATED_HITS" && rm -f "${ANNOTATED_HITS}.bak"
            echo "Warning: Species annotation failed - using unannotated results"
        fi
    else
        cp "$TOP_HITS" "$ANNOTATED_HITS"
        sed -i.bak '1s/$/ species/' "$ANNOTATED_HITS" && rm -f "${ANNOTATED_HITS}.bak"
    fi
   
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
        awk -F'\t' '{printf "%-15s %-15s %-6s %-6s %-8s %-7s %-6s %-6s %-6s %-6s %-10s %-10s %s\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' "$ANNOTATED_HITS"
    ) > "$REPORT_FILE"
   
    echo " - Results saved to $SAMPLE_DIR"
done

echo -e "\nProcessing complete!"
echo "Output directories created:"
find "$OUTPUT_DIR" -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
    echo "- $dir"
done

exit 0
