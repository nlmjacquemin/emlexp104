#!/bin/bash

# Default parameters
DEFAULT_LOG_FOLDER=~/jobs
DEFAULT_WORKING_DIRECTORY=/data
DEFAULT_INPUT_FOLDER="raw/reads"
DEFAULT_OUTPUT_FOLDER="analysis/fastqc/raw"
DEFAULT_INPUT_EXTENSION=".fastq"  # You can set it to .fq or other as necessary

# Show usage function
show_usage() {
    echo "Usage: $0"
    echo "Options:"
    echo "  -w, --working_directory  Specify the path of the working directory (default: $DEFAULT_WORKING_DIRECTORY)"
    echo "  -i, --input_folder       Specify the path of the input read folder (default: $DEFAULT_INPUT_FOLDER)"
    echo "  -x, --input_extension    Specify the input reads extension (default: $DEFAULT_INPUT_EXTENSION)"
    echo "  -o, --output_folder      Specify the path of the output folder (default: $DEFAULT_OUTPUT_FOLDER)"
    echo "  -l, --log_folder         Specify the log folder (default: $DEFAULT_LOG_FOLDER)"
    echo "  --                       Separator to pass following options to the script/command being executed"
}

# Default variables
working_directory="$DEFAULT_WORKING_DIRECTORY"
input_folder="$DEFAULT_INPUT_FOLDER"
output_folder="$DEFAULT_OUTPUT_FOLDER"
input_extension="$DEFAULT_INPUT_EXTENSION"
log_folder="$DEFAULT_LOG_FOLDER"

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
        -x|--input_extension)
            input_extension="$2"
            shift 2
            ;;
        -o|--output_folder)
            output_folder="$2"
            shift 2
            ;;
        -l|--log_folder)
            log_folder="$2"
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

cmd_options=("$@")

# Function to print log messages with timestamps
log_message() {
    local message=$@
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "[${timestamp}] ${message}"
}

# Check if fastqc is installed
if ! command -v fastqc &> /dev/null; then
    echo "fastqc could not be found. Please install it and ensure it's in your PATH."
    exit 1
fi

# Log start of the process
log_message "Starting FastQC analysis in $working_directory"

# Move to the working directory
cd "$working_directory" || exit

# Create output folder
mkdir -p "$output_folder"

# Process all the input files at once
log_message "Processing all FASTQ files in $input_folder"

fastqc \
    "$input_folder"/*"${input_extension}" \
    -o "$output_folder" \
    --threads "$SLURM_CPUS_PER_TASK" \
    "${cmd_options[@]}"

# Run MultiQC after FastQC is done
log_message "Running MultiQC..."
multiqc \
    --outdir "$output_folder" \
    "$output_folder"  \
    "${cmd_options[@]}"

log_message "FastQC analysis completed!"
log_message "QC can be found in $output_folder"
