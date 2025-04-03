####--Library####

library(optparse)
library(Biostrings)
library(dada2)
library(stringr)
library(rjson) #Json format

sessionInfo()

####--Getting today date####

date<-Sys.Date()
time<-format(Sys.time(), "%H-%M-%S")

####--Detecting arguments####

args_list = list(
  make_option(
    c("--debug"),
    type = "logical",
    default = FALSE, 
    help = "Enable debug mode."
  ),
  make_option(
    c("--config"),
    type = "character",
    default = "config.json", 
    help = "Path to the config file either relative to working directory or absolute."
  ),
  make_option(
    c("--db_folder"),
    type = "character",
    default = "db",
    help = "Path to folder containing databases."
  )
  # make_option(
  #   c("--pooling"),
  #   type = "character",
  #   default = FALSE, 
  #   help = "If pooling = TRUE, the algorithm will pool together all samples prior to sample inference. If pooling = FALSE, sample inference is performed on each sample individually. If pooling = pseudo, the algorithm will perform pseudo-pooling between individually processed samples."
  # )
) 

arg_parser = OptionParser(option_list=args_list)
args = parse_args(arg_parser)

print(args)

####--Setting working directory####

if (args$debug){
  
  setwd("C:/Users/nljacque/Desktop/Share_VM_Ubuntu/container/1_Project/0_Commun/eml_scitas/tests/result/pacbio/dada")
  
} 

#By default, the working directory is where the script have been called by Rscript 

print(getwd())

####--Loading configuration####

if (file.exists(args$config)){
  
  config <- fromJSON(file = args$config)
  
} else {
  
  config<-data.frame()
  
}

print("Configuration:")

print(config)

####--Preparing input and output####

# For abundance to taxonomy (dechimerization, taxonomic assignment)

## Setting Folder

if (is.null(config$output_folder)){
  
  input_folder <- "dada/result"
  
} else {
  
  input_folder  <- config$output_folder
  
}

if (is.null(config$output_folder)){
  
  output_folder <- "picrust/raw"
  
} else {
  
  output_folder <- config$output_folder
  
}

input_intermediary<-file.path(
  input_folder,
  "intermediary"
)

## Create if not existing

if (!dir.exists(output_folder)){
  
  dir.create(
    output_folder,
    recursive=T
    )
  
}

####--Transforming data####

# Abundance table

## Loading data

abd<-readRDS(
  file.path(
    input_intermediary,
    "abd_dechimerized.rds"
    )
  )

## Saving data 

ASV_ids <- paste0(
  "ASV",
  seq(ncol(abd))
  )

abd_wide<-as.data.frame(
  cbind(
    data.frame( ASV_id = ASV_ids),
    t(abd),
    row.names = NULL
  )
)

print(head(abd_wide))

write.table(
  x = abd_wide,
  file = file.path(
    input_intermediary,
    "abd_uncorrected.wide.tsv"
  ),
  sep="\t",
  row.names = F
)

file.copy(
  file.path(
    input_intermediary,
    "abd_uncorrected.wide.tsv"
  ),
  file.path(
    output_folder,
    "abd_uncorrected.wide.tsv"
  )
)

# Taxonomic table

tax_files<-list.files(
  path = file.path(
    input_intermediary
  ),
  pattern = "tax_.*\\.rds"
)

print(tax_files)

for (tax_file in tax_files){

  # Loading data

  tax<-readRDS(
    file.path(
      input_intermediary,
      tax_file
      )
    )

  tax_wide<-as.data.frame(
    cbind(
      data.frame( ASV_id = ASV_ids),
      tax,
      row.names = NULL
    )
  )

  print(head(tax))
 
  write.table(
    x = tax_wide,
    file = file.path(
      input_intermediary,
      str_replace(tax_file,"\\.rds","\\.tsv")
    ),
    sep="\t",
    row.names = F
  )

}
# Sequence Fasta

sequences <- DNAStringSet(
  colnames(abd)
) # Create a DNAStringSet from the ASVs

names(sequences) <- ASV_ids

writeXStringSet(
  sequences,
  file.path(
    output_folder,
    "ASV.fna"
  )
)

seq_wide<-as.data.frame(
  cbind(
    data.frame( ASV_id = ASV_ids),
    sequences,
    row.names = NULL
  )
)

write.table(
  x = seq_wide,
  file = file.path(
    input_intermediary,
    "seq.tsv"
  ),
  sep="\t",
  row.names = F
)

