#!/bin/bash

# Default values
DEFAULT_LOG_FOLDER=~/jobs
DEFAULT_WORKING_DIRECTORY=/data
DEFAULT_SCRIPT_FOLDER=~/scripts
DEFAULT_SCRIPT_EXTENSION=".R"

# Show usage information
show_usage() {
    echo "Usage: $0"
    echo "Options:"
    echo "  -w, --working_directory     Specify the path of the working directory (default: $DEFAULT_WORKING_DIRECTORY)"
    echo "  -s, --script_folder         Specify the path of the script folder (default: $DEFAULT_SCRIPT_FOLDER)"
    echo "  -r, --r_script              Specify the R script to run"
    echo "  -x, --script_extension      Specify the extension (default: $DEFAULT_SCRIPT_EXTENSION)"
    echo "  -l, --log_folder            Specify the log folder (default: $DEFAULT_LOG_FOLDER)"
    echo "  --                          Separator to pass following options to the R script"
    echo "  [r_script_options]          Options for the R script"
}

# Default variables
working_directory="$DEFAULT_WORKING_DIRECTORY"
script_folder="$DEFAULT_SCRIPT_FOLDER"
r_script=""
script_extension=$DEFAULT_SCRIPT_EXTENSION
log_folder="$DEFAULT_LOG_FOLDER"

# Parse options
while [[ "$1" != "" ]]; do
    case "$1" in
        -w|--working_directory)
            working_directory="$2"
            shift 2
            ;;
        -s|--script_folder)
            script_folder="$2"
            shift 2
            ;;
        -r|--r_script)
            r_script="$2"
            shift 2
            ;;
        -x|--script_extension)
            script_extension="$2"
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

# Remaining arguments are for the R script
r_script_options=("$@")

# Function to print log messages with timestamp
log_message() {
    local message=$@
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "[${timestamp}] ${message}"
}

# Change to the working directory
log_message "Changing to working directory: $working_directory"
cd "$working_directory"

# Path to the R script
path_to_script="$script_folder/${r_script}${script_extension}"

# Check if the R script exists
if [[ ! -f "$path_to_script" ]]; then
    log_message "R script not found: $path_to_script"
    exit 1
fi

# Run the R script
log_message "Running R script: $path_to_script"
Rscript "$path_to_script" "${r_script_options[@]}"
log_message "R script execution completed"

# Copy the R script to the log folder
log_message "Copying R script to log folder: $log_folder"
mkdir -p "$log_folder"
cp "$path_to_script" "$log_folder"
log_message "Copying completed"
