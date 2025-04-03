#!/bin/bash

#SBATCH -vv
#SBATCH --nodes 1
#SBATCH --mail-user nicolas.jacquemin@epfl.ch
#SBATCH --ntasks 1
#SBATCH --chdir /home/nljacque/jobs/tmp
#SBATCH --output slurm-launcher-%j.out

SCRIPT_VERSION=2.1.0

DEFAULT_SCRIPT_FOLDER=~/scripts
DEFAULT_SCRIPT_EXTENSION=".sh"
DEFAULT_LOG_FOLDER=~/jobs
DEFAULT_WORKING_LOG_FOLDER=scripts

show_usage() {
    echo "Usage: $0 [-s script_name] [-f script_folder] [-l log_folder] -- [script_options]"
    echo ""
    echo "Options:"
    echo "	-s, --script_name           Specify the script name to launch (without extension)"
    echo "	-x, --script_extension       Specify extension of the script (default: $DEFAULT_SCRIPT_EXTENSION)"
    echo "	-f, --script_folder          Specify the path of the folder containing the script (default: $DEFAULT_SCRIPT_FOLDER)"
    echo "	-l, --log_folder             Specify the log folder (default: $DEFAULT_LOG_FOLDER)"
    echo "  -w, --working_directory      Specify the path of the working directory"
    echo "  -k, --working_log_folder     If used, create a log folder in the working directory at the specified path (default: $DEFAULT_WORKING_LOG_FOLDER)"
    echo "	--                           Separator to pass following options to the script being launched"
    echo "	[script_options]             Options for the script being launched"
    echo "  -h, --help                   Show this help message"
}

# Default values
script_name=""
script_extension="$DEFAULT_SCRIPT_EXTENSION"
script_folder="$DEFAULT_SCRIPT_FOLDER"
log_folder="$DEFAULT_LOG_FOLDER"
working_log_folder="$DEFAULT_WORKING_LOG_FOLDER"

# Parse options
while [[ "$1" != "" ]]; do
    case "$1" in
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
        -w|--working_directory)
            working_directory="$2"
            shift 2
            ;;
        -k|--working_log_folder)
            if [[ -n "$2" && "$2" != -* ]]; then
                working_log_folder="$2"
                shift 2
            else
                shift 1
            fi
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

# Remaining arguments are for the script to run
script_options=("$@")

# Validate required parameters
if [[ -z "$script_name" ]]; then
    echo "Error: Script name must be specified with -s."
    show_help
    exit 1
fi

# Function to print log messages with timestamp
log_message() {
    local message=$@
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    echo -e "[${timestamp}] ${message}"
}

# Set unique log directory for each job
if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
    job_log_folder="$log_folder/$SLURM_ARRAY_JOB_ID/$SLURM_ARRAY_TASK_ID"
    working_log_folder="$working_log_folder/$SLURM_ARRAY_JOB_ID"  # Append array job ID to working log folder
else
    job_log_folder="$log_folder/$SLURM_JOB_ID"
fi
mkdir -p "$job_log_folder"

# Print job info for logs
log_message "Setting batch directives..."
echo -e "Job Name: $SLURM_JOB_NAME"
echo -e "CPUs per Task: $SLURM_CPUS_PER_TASK"
echo -e "Total Tasks: $SLURM_NTASKS"
echo -e "Job ID: $SLURM_JOBID"
if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
    echo -e "Array Job ID: $SLURM_ARRAY_JOB_ID"
    echo -e "Array Task ID: $SLURM_ARRAY_TASK_ID"
fi
log_message "Setting batch directives finished!"

# Setting logs
log_message "Setting logs..."
path_to_script="$script_folder/${script_name}${script_extension}"
cp "$(realpath $0)" "$job_log_folder/launcher.sh"
cp "$path_to_script" "$job_log_folder"
log_message "Setting logs finished!"

# Launching the job with Apptainer
log_message "Launching the job with Apptainer."
srun -J "$script_name" \
    --output="$job_log_folder/slurm-%j.out" \
    bash "$path_to_script" --log_folder "$job_log_folder" --working_directory "$working_directory" "${script_options[@]}"
log_message "The job is finished!"

# Moving log files
log_message "Moving log files..."
mv /home/nljacque/jobs/tmp/slurm*${SLURM_JOB_ID}* "$job_log_folder"
if [[ -n "$working_log_folder" ]]; then
    mkdir -p "$working_directory/$working_log_folder"
    cp -r "$job_log_folder" "$working_directory/$working_log_folder"
    
    if [[ -n "$SLURM_ARRAY_TASK_ID" ]]; then
        echo -e "Log files can be found in: $job_log_folder and $working_directory/$working_log_folder/${SLURM_ARRAY_TASK_ID}"
    else
        echo -e "Log files can be found in: $job_log_folder and $working_directory/$working_log_folder/${SLURM_JOB_ID}"
    fi
else
    echo -e "Log files can be found in: $job_log_folder"
fi

log_message "Moving log files finished!"
