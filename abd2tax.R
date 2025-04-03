#### --Library####

library(optparse)
library(dada2)
library(stringr)
library(rjson) # Json format

sessionInfo()

#### --Getting today date####

date <- Sys.Date()
time <- format(Sys.time(), "%H-%M-%S")

#### --Detecting arguments####

args_list <- list(
  make_option(
    c("--debug"),
    type = "logical",
    default = FALSE,
    help = "Enable debug mode."
  ),
  make_option(
    c("--threads"),
    type = "integer",
    default = 1,
    help = "Number of threads to use."
  ),
  make_option(
    c("--dowloading_db"),
    type = "logical",
    default = FALSE,
    help = "IF database should be downloaded."
  ),
  make_option(
    c("--db_folder"),
    type = "character",
    default = "db",
    help = "Path to folder containing databases."
  ),
  make_option(
    c("--config"),
    type = "character",
    default = "config.json",
    help = "Path to the config file either relative to working directory or absolute."
  )
  # make_option(
  #   c("--pooling"),
  #   type = "character",
  #   default = FALSE,
  #   help = "If pooling = TRUE, the algorithm will pool together all samples prior to sample inference. If pooling = FALSE, sample inference is performed on each sample individually. If pooling = pseudo, the algorithm will perform pseudo-pooling between individually processed samples."
  # )
)

arg_parser <- OptionParser(option_list = args_list)
args <- parse_args(arg_parser)

print(args)

#### --Setting working directory####

if (args$debug) {
  setwd("C:/Users/nljacque/Desktop/Share_VM_Ubuntu/container/1_Project/0_Commun/eml_scitas/tests/result/pacbio/dada")
}

# By default, the working directory is where the script have been called by Rscript

print(getwd())

#### --Function####

load_abd <- function(
    input_folder = "result/intermediary",
    input_file = "abd_undechimerized.rds") {
  input_file_path <- file.path(
    input_folder,
    input_file
  )

  if (file.exists(input_file_path)) {
    print(paste0("Loading and merging:", input_file_path))
    abd <- readRDS(input_file_path)
  } else {
    stop(paste0("File not found: ", input_file_path, ". Stopping execution."))
  }

  return(abd)
}

#### --Loading configuration####

if (file.exists(args$config)) {
  config <- fromJSON(file = args$config)
} else {
  config <- NULL
}

print("Configuration:")

print(config)

#### --Preparing input and output####

# For abundance to taxonomy (dechimerization, taxonomic assignment)

## Setting Folder

if (is.null(config$output_folder)) {
  input_folder <- "result"
} else {
  input_folder <- config$output_folder
}

input_intermediary <- file.path(
  input_folder,
  "intermediary"
)


if (is.null(config$output_folder)) {
  output_folder <- "result"
} else {
  output_folder <- config$output_folder
}

output_intermediary <- file.path(
  output_folder,
  "intermediary"
)

## Create if not existing

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

output_intermediary <- file.path(
  output_folder,
  "intermediary"
)

if (!dir.exists(output_intermediary)) {
  dir.create(output_intermediary)
}

#### --Loading abundance files####

if (is.null(config$abd_default_file)) {
  abd_default_file <- "abd_undechimerized.rds"
} else {
  abd_default_file <- config$abd_default_file
}

abd <- load_abd(
  input_folder = input_intermediary,
  input_file = abd_default_file
)

rownames(abd) <- str_remove(rownames(abd), ".gz|.fastq")

#### --Removing chimera####

print("Starting removing chimeras")

abd <- removeBimeraDenovo(
  abd,
  method = "consensus", # "consensus"(default),"pooled" or "per-sample"
  multithread = args$threads
)

print("Removing chimeras done!")

write.table(
  x = t(abd),
  file = file.path(
    output_intermediary,
    "abd_dechimerized.tsv"
  ),
  sep = "\t",
  row.names = T
)

saveRDS(
  abd,
  file.path(
    output_intermediary,
    "abd_dechimerized.rds"
  )
)

#### --Taxonomic assignment####

# Set up databases folder

if (!dir.exists(args$db_folder)) {
  stop(paste0("Error: No ", args$db_folder, " folder found at ", getwd()))
}

if (is.null(config$databases)) {
  dbs <- c("silva")
} else {
  dbs <- config$databases
}

# Download database

if ("silva" %in% dbs) {
  db <- "silva"

  one_db_folder <- file.path(
    args$db_folder,
    db
  )

  one_db_folder_format <- file.path(
    one_db_folder,
    "dada_format"
  )

  if (args$dowloading_db) {
    dir.create(one_db_folder)

    print(paste0("Downloading silva database at:", one_db_folder))

    download.file(
      "https://www.arb-silva.de/fileadmin/silva_databases/release_138/Exports/SILVA_138_SSURef_tax_silva.fasta.gz",
      file.path(one_db_folder, "SILVA_138_SSURef_tax_silva.fasta.gz")
    )

    dir.create(one_db_folder_format)

    print(paste0("Downloading silva database in dada_format at:", one_db_folder_format))

    # For short reads (assignTaxonomy + addSpecies)

    ## To be used by assignTaxonomy

    download.file(
      "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_train_set.fa.gz",
      file.path(one_db_folder_format, "SILVA_138_SSURef_tax_silva.fasta.gz")
    )

    ## To be used by addSpecies

    download.file(
      "https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz",
      file.path(one_db_folder_format, "silva_species_assignment_v138.1.fa.gz")
    )

    # For long reads (assignTaxonomy)

    ## To be used by assignTaxonomy

    download.file(
      "https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz",
      file.path(one_db_folder_format, "silva_nr99_v138.1_wSpecies_train_set.fa.gz")
    )
  }

  # Taxonomic assignment

  print("Starting taxonomic assignment")

  tax <- assignTaxonomy(
    seqs = abd,
    refFasta = file.path(
      one_db_folder_format,
      "silva_nr99_v138.1_wSpecies_train_set.fa.gz"
    ),
    multithread = args$threads
  )

  print("Taxonomic assignment done!")

  write.table(
    x = tax,
    file = file.path(
      output_intermediary,
      paste0(
        "tax_",
        db,
        ".tsv"
      )
    ),
    sep = "\t",
    row.names = T
  )

  saveRDS(
    tax,
    file.path(
      output_intermediary,
      paste0(
        "tax_",
        db,
        ".rds"
      )
    )
  )
}

#### --Intermediary saving####

print("Saving R workspace...")

save.image(
  file = file.path(
    output_intermediary,
    "abd2tax.RData"
  )
)

print("Saving R workspace done!")
