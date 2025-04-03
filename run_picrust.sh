#!/bin/bash

# Default values
DEFAULT_LOG_FOLDER=~/jobs
DEFAULT_WORKING_DIRECTORY=/data
DEFAULT_INPUT_FOLDER="raw"
DEFAULT_OUTPUT_FOLDER="result"
DEFAULT_SEQUENCE_FILE="ASV.fna"
DEFAULT_ABUNDANCE_FILE="abd_uncorrected.wide.tsv"
DEFAULT_CPUS=${SLURM_CPUS_PER_TASK:-4}

# Show usage information
show_usage() {
    echo "Usage: $0"
    echo "Options:"
    echo "  -w, --working_directory    Specify the path of the working directory (default: $DEFAULT_WORKING_DIRECTORY)"
    echo "  -i, --input_folder         Specify the input folder (default: $DEFAULT_INPUT_FOLDER)"
    echo "  -o, --output_folder        Specify the output folder (default: $DEFAULT_OUTPUT_FOLDER)"
    echo "  -s, --sequence_file        Specify the sequence file (default: $DEFAULT_SEQUENCE_FILE)"
    echo "  -a, --abundance_file       Specify the abundance file (default: $DEFAULT_ABUNDANCE_FILE)"
    echo "  -l, --log_folder           Specify the log folder (default: $DEFAULT_LOG_FOLDER)"
    echo "  -p, --cpus                 Specify the number of CPUs to use (default: $DEFAULT_CPUS)"
    echo "  --                         Separator to pass following options to picrust2_pipeline.py"
}

# Default variables
log_folder="$DEFAULT_LOG_FOLDER"
working_directory="$DEFAULT_WORKING_DIRECTORY"
input_folder="$DEFAULT_INPUT_FOLDER"
output_folder="$DEFAULT_OUTPUT_FOLDER"
sequence_file="$DEFAULT_SEQUENCE_FILE"
abundance_file="$DEFAULT_ABUNDANCE_FILE"
cpus="$DEFAULT_CPUS"

# Parse options
while [[ "$1" != "" ]]; do
    case "$1" in
        -w|--working_directory)
            working_directory="$2"
            shift 2
            ;;
        -i|--input_folder)
            input_folder="$2"
            shift 2
            ;;
        -o|--output_folder)
            output_folder="$2"
            shift 2
            ;;
        -s|--sequence_file)
            sequence_file="$2"
            shift 2
            ;;
        -a|--abundance_file)
            abundance_file="$2"
            shift 2
            ;;
        -l|--log_folder)
            log_folder="$2"
            shift 2
            ;;
        -p|--cpus)
            cpus="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
    esac
done

# Remaining arguments are passed to picrust2_pipeline.py
cmd_options=("$@")

# Log setup
day=$(date +%F)
hour=$(date +%T)
log_file="$log_folder/log_picrust_${day}_$hour.txt"

# Print log message function
log_message() {
    local message=$@
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "[${timestamp}] ${message}"
}

log_message "Starting PICRUSt2 pipeline"

# Change to working directory
cd "$working_directory"

# Run PICRUSt2 pipeline
picrust2_pipeline.py \
    --stratified --per_sequence_contrib \
    -s "$input_folder/$sequence_file" \
    -i "$input_folder/$abundance_file" \
    -o "$output_folder" \
    -p "$cpus" \
    --verbose "${cmd_options[@]}" \
    2>&1 | tee "$log_file"

# Add descriptions for EC, KO, and METACYC
log_message "Adding descriptions to output files"

add_descriptions.py -i "$output_folder/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" -m EC -o "$output_folder/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz"
add_descriptions.py -i "$output_folder/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" -m KO -o "$output_folder/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz"
add_descriptions.py -i "$output_folder/pathways_out/path_abun_unstrat.tsv.gz" -m METACYC -o "$output_folder/pathways_out/path_abun_unstrat_descrip.tsv.gz"

log_message "PICRUSt2 pipeline finished"
