#!/bin/bash

SCRIPT_VERSION=2.0.0

DEFAULT_IMAGE_FOLDER=~/images
DEFAULT_IMAGE_EXTENSION=".sif"
DEFAULT_SCRIPT_FOLDER=~/scripts
DEFAULT_SCRIPT_EXTENSION=".sh"
DEFAULT_LOG_FOLDER=~/jobs
DEFAULT_CONTAINER_WORKING_DIRECTORY=/data

show_usage() {
    echo "Usage: $0 [-w working_directory] [-i image_name] [-s script_name] [-l log_folder] [--micromamba_env environment_name] [apptainer_options] -- [script_options]"
    echo ""
    echo "This script sets up and execute an Apptainer (formerly Singularity) container environment, executing a specified script within the container." 
    echo "It allows for configuration of the working directory, image, script, and log folder, as well as additional Apptainer options."
    echo ""
    echo "Options:"
    echo "  -w, --working_directory  Specify the path in the host of the working directory"
    echo "  -i, --image_name         Specify the image name (without extension)"
    echo "  -y, --image_extension    Specify the extension of the image (default: $DEFAULT_IMAGE_EXTENSION)"
    echo "  -d, --image_folder       Specify the path of the folder containing the images (default: $DEFAULT_IMAGE_FOLDER)"
    echo "  -s, --script_name        Specify the script name to execute (without extension)"
    echo "  -x, --script_extension   Specify the extension of the script (default: $DEFAULT_SCRIPT_EXTENSION)"
    echo "  -f, --script_folder      Specify the path of the folder containing the script (default: $DEFAULT_SCRIPT_FOLDER)"
    echo "  -l, --log_folder         Specify the log folder (default: $DEFAULT_LOG_FOLDER)"
    echo "  --micromamba_env         Specify the Micromamba environment name"
    echo "  [apptainer_options]      Options for Apptainer exec command"
    echo "  --                       Separator to pass following options to the script being execute"
    echo "  [script_options]         Options for the script being execute (optional)"
    echo "  -h, --help               Show this help message"
    echo ""
    echo "Note: A binding is by default set up for the host working directory inside the container as $DEFAULT_CONTAINER_WORKING_DIRECTORY"
}

# Fetch Apptainer exec options dynamically
fetch_apptainer_options() {
    apptainer exec --help | grep -oP '(?<=\s)--\w+(-\w+)?' | uniq
}

# Validate Apptainer option
is_apptainer_option() {
    local opt="$1"
    # Fetch Apptainer options once at the start
    apptainer_options_list=($(fetch_apptainer_options))
    for apptainer_opt in "${apptainer_options_list[@]}"; do
        if [[ "$opt" == "$apptainer_opt" ]]; then
            return 0
        fi
    done
    return 1
}

# Default values
path_to_host_working_directory=""
path_to_container_working_directory="$DEFAULT_CONTAINER_WORKING_DIRECTORY"
image_name=""
image_extension="$DEFAULT_IMAGE_EXTENSION"
image_folder="$DEFAULT_IMAGE_FOLDER"
script_name=""
script_extension="$DEFAULT_SCRIPT_EXTENSION"
script_folder="$DEFAULT_SCRIPT_FOLDER"
log_folder="$DEFAULT_LOG_FOLDER"
micromamba_env=""
apptainer_options=()

# Parse options
while [[ "$1" != "" ]]; do
    case "$1" in
        -w|--working_directory)
            path_to_host_working_directory="$2"
            shift 2
            ;;
        -i|--image_name)
            image_name="$2"
            shift 2
            ;;
        -y|--image_extension)
            image_extension="$2"
            shift 2
            ;;
        -d|--image_folder)
            image_folder="$2"
            shift 2
            ;;
        -s|--script_name)
            script_name="$2"
            shift 2
            ;;
        -x|--script_extension)
            script_extension="$2"
            shift 2
            ;;
        -f|--script_folder)
            script_folder="$2"
            shift 2
            ;;
        -l|--log_folder)
            log_folder="$2"
            shift 2
            ;;
        --micromamba_env)
            micromamba_env="$2"
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
            if is_apptainer_option "$1"; then
                apptainer_options+=("$1")
                shift
                # If the option takes a value, include it as well
                if [[ "$1" != "" && ! "$1" =~ ^- ]]; then
                    apptainer_options+=("$1")
                    shift
                fi
            else
                echo "Unknown option: $1"
                show_usage
                exit 1
            fi
            ;;
    esac
done

# Remaining arguments are for the script to run
script_options=("$@")

# Validate required parameters
if [[ -z "$path_to_host_working_directory" || -z "$image_name" || -z "$script_name" ]]; then
    echo "Error: Missing mandatory options."
    show_help
    exit 1
fi

# Function to print log messages with timestamp
log_message() {
    local message=$@
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "[${timestamp}] ${message}"
}

# Display information for logs
log_message "Setting apptainer container..."
# User
USERNAME=$(whoami)
echo "USERNAME:" $USERNAME
# Working directory
echo "path_to_host_working_directory:" $path_to_host_working_directory
# Image
path_to_the_image="$image_folder/${image_name}${image_extension}"
echo "image:" $path_to_the_image
# Script
path_to_script_host="$script_folder/${script_name}${script_extension}"
path_to_script_container=$path_to_script_host
echo "script:" $path_to_script_container
# Saving the script in log folder
echo "path_to_log_directory:" $log_folder
cp $path_to_script_host $log_folder
log_message "Setting apptainer container finished!"

#Executing script within apptainer container
log_message "Running apptainer container and executing the script..."
apptainer_cmd=(
    --bind "$path_to_host_working_directory":"$path_to_container_working_directory"
    "${apptainer_options[@]}"
    "$path_to_the_image"
)
script_cmd=(
    bash "$path_to_script_container"
    --working_directory "$path_to_container_working_directory"
    --log_folder "$log_folder"
    "${script_options[@]}"
)
if [[ -n "$micromamba_env" ]]; then
    log_message "\nActivating micromamba environment..."
    echo "micromamba_env:" "$micromamba_env"
    apptainer exec \
        "${apptainer_cmd[@]}" \
        micromamba run -n "$micromamba_env" \
        "${script_cmd[@]}"
    log_message "Closing micromamba environment."
else
    apptainer exec \
        "${apptainer_cmd[@]}" \
        "${script_cmd[@]}"
fi
log_message "Running apptainer container and executing the script done!"
