#### --Library####

library(optparse)
library(dada2)
library(stringr)
library(rjson) # Json format
library(dplyr) # %>% mutate etc..
library(reshape2) # melt
library(scales) # scientific_format
library(ggplot2)
# library(catecolors) # glasbey()

#### --Getting date and time####

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
    c("--big_data"),
    type = "logical",
    default = FALSE,
    help = "If the data should consider as big data or not. If you want to do pooling for denoising let the value by default ('FALSE')."
  ),
  make_option(
    c("--threads"),
    type = "integer",
    default = 1,
    help = "Number of threads to use."
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

args_table <- as.data.frame(
  args,
  row.names = c("Arguments")
)

print(
  args_table
)

#### --Setting working directory####

if (args$debug) {
  setwd("C:/Users/nljacque/Desktop/Share_VM_Ubuntu/container/1_Project/0_Commun/eml_scitas/tests/result/pacbio/dada")
}

# By default, the working directory is where the script have been called by Rscript

print(getwd())

#### --Function####

# Custom label

custom_label <- function(x) {
  parse(text = gsub("e", " %*% 10^", scientific_format()(x)))
}

# Plotting read length distribution

plot_read_hist <- function(
    files,
    step_name,
    ...) {
  kwargs <- list(...)

  if (is.null(kwargs$output_folder)) {
    kwargs$output_folder <- getwd()
  }

  dl <- lapply(files, function(x) getSequences(x))

  names(dl) <- c(samples)

  df <- melt(dl)

  colnames(df) <- c("sequences", "samples")

  df <- df %>%
    mutate(df, length = nchar(sequences))

  p <- ggplot(
    df,
    aes(
      x = length
    )
  ) +
    geom_histogram(
      aes(
        fill = samples
      )
    ) +
    stat_bin(
      geom = "text",
      aes(
        label = ifelse(
          after_stat(count) != 0,
          after_stat(count),
          ""
        )
      ),
      vjust = -0.5,
      color = "black",
      size = 3,
      fontface = "bold"
    ) +
    xlab("Read Length") +
    ylab("Read Count") +
    # scale_fill_manual(
    #  values = glasbey()[-1]
    #  ) +
    scale_y_continuous(
      labels = custom_label
    ) +
    theme(legend.position = "none")
  # coord_trans(y = "pseudo_log") #/!\ Not correct because transformation of fill is non linear
  # coord_cartesian(ylim = c(0, 10))

  ggsave(
    filename = file.path(
      output_folder,
      paste0(
        "Read_length_distribution_after_",
        step_name,
        ".png"
      )
    ),
    plot = p, # By default plot = last_plot()
    scale = 1,
    width = 20,
    height = 20,
    units = "cm"
  )

  return(p)
}

# Saving read amounts tracked

save_track <- function(
    df,
    step_name,
    ...) {
  kwargs <- list(...)

  if (is.null(kwargs$output_folder)) {
    kwargs$output_folder <- getwd()
  }

  print("Saving read tracking...")

  file <- file.path(
    kwargs$output_folder,
    "read_track.tsv"
  )

  df <- as.data.frame(df)

  colnames(df) <- c("raw", step_name)

  if (file.exists(file)) {
    print(
      paste0(
        "Append pre-existing read tracking file at: ",
        file
      )
    )

    pre_df <- read.table(
      file = file,
      header = TRUE,
      sep = "\t"
    )

    if (!all(rownames(df) == pre_df$samples)) {
      stop("Error: Pre-existing read_track.tsv file that does not have the same samples or samples ordered in the same way.")
    }

    df <- cbind(
      pre_df,
      df[-1]
    )
  } else {
    df <- cbind(
      data.frame(
        samples = rownames(df)
      ),
      df
    )
  }

  write.table(
    x = df,
    file = file,
    sep = "\t",
    row.names = FALSE
  )

  print(
    paste0(
      "Read tracking saved at: ",
      file
    )
  )
}

# Sanity check

sanity_check <- function(
    files,
    step_name,
    track,
    output_folder,
    ...) {
  save_track(
    df = as.data.frame(track),
    step_name = step_name,
    output_folder = output_folder
  )

  plot_read_hist(
    files = files,
    step_name = step_name,
    output_folder = output_folder
  )
}

# Saving function

saving_data <- function(
    data,
    data_name,
    output_folder,
    seq_col,
    track = NULL) {
  if (is.null(track)) {
    track <- data.frame(
      nothing = rep(
        NA,
        length(data)
      ),
      raw = rep(
        NA,
        length(data)
      )
    )
  }

  track <- data.frame(
    row.names = samples,
    raw = track[, 2],
    data.frame(
      dereplicated = sapply(data, function(x) sum(x[[seq_col]]))
    )
  )

  save_track(
    df = as.data.frame(track),
    step_name = data_name,
    output_folder = output_folder
  )

  # data_name = deparse(substitute(data))

  file_out <- file.path(
    output_folder,
    paste0(
      data_name,
      ".rds"
    )
  )

  saveRDS(
    data,
    file = file_out
  )

  print(
    paste0(
      "Reads ",
      data_name,
      " saved at:",
      file_out
    )
  )
}

learn_errors_step <- function(
    derep,
    output_folder,
    output_intermediary,
    threads) {
  print("Starting learning error rates")

  set.seed(100)

  err <- learnErrors(
    derep,
    BAND_SIZE = 32,
    randomize = TRUE,
    multithread = threads,
    errorEstimationFunction = PacBioErrfun # loessErrfun by default
    # nbases = 1e+08 #by default nbases = 1e+08
  )
  # /!\ Should use #PacBioErrfun /!\

  print("Learning error rates done!")

  plotErrors(
    err,
    nominalQ = T
  )

  # Reminder for nominalQ:
  # Q = -10*log10(P) <=> P = 10^(-Q/10)
  # Q : phred/consensus quality score
  # P : probabilities of base-calling error
  # The probabilities of base being wrongly assign to ONE other base is 1/3 of Q !
  # The base A can be wrongly assigned between the bases C,T,G each wrong calling sharing 1/3 of the total error rate.

  ggsave(
    paste0(output_folder, "/Error_rates.png"),
    width = 20,
    height = 20,
    units = "cm"
  )

  data <- err

  data_name <- "errors"

  file_out <- file.path(
    output_intermediary,
    paste0(
      data_name,
      ".RData"
    )
  )

  save(
    data,
    file = file_out
  )

  print(
    paste0(
      "Error rates saved at:",
      file_out
    )
  )

  # saving_data(
  #   data = err,
  #   data_name = "errors",
  #   output_folder = output_intermediary,
  #   seq_col = "uniques",
  #   track = track
  # )

  return(err)
}

#### --Loading configuration file####

if (file.exists(args$config)) {
  config <- fromJSON(file = args$config)
} else {
  config <- NULL
}

print("Configuration:")

print(config)

# Pooling options

if (is.null(config$pooling)) {
  pooling <- FALSE
} else {
  pooling <- config$pooling
}

# Warning

if (pooling %in% c("pseudo", TRUE) & args$big_data == TRUE) {
  cat(
    paste0(
      "'Pseudo' or 'TRUE' pooling would not be effective when big_data is 'TRUE' as samples are processed one by one.",
      "\nSet big_data to 'FALSE' if you want to do denoising with 'pseudo' or 'TRUE' pooling."
    )
  )
}

#### --Preparing input and output####

# For read to dereplicates (removing primers, filtering, dereplication)

## Setting Folder

if (is.null(config$input_folder)) {
  input_folder <- "raw/reads"
} else {
  input_folder <- config$input_folder
}

if (is.null(config$output_folder)) {
  output_folder <- "result"
} else {
  output_folder <- config$output_folder
}

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

#### --Saving arguments inputed####

write(
  toJSON(config),
  file = file.path(output_intermediary, "config.json")
)

write.csv(
  x = args_table,
  file = file.path(output_intermediary, "args.tsv"),
  sep = "\t"
)

#### --Loading read files####

input_files <- file.path(
  input_folder,
  list.files(input_folder)[
    grepl(
      ".fastq|.fastq.gz",
      list.files(input_folder)
    )
  ]
)

samples <- str_remove_all(basename(input_files), ".fastq|.hifi_reads|demultiplex.|.gz")

print(paste("Loading", length(input_files), "sample fastq files."))

print("Samples:")

cat(paste(collapse = "\n", samples))

#### --Removing Primers####

# Getting primers sequences

print("Note: The forward and reverse primer sequences depends on the kit used.")

# Forward primer

if (is.null(config$primers$fwd)) {
  print("No forward primer sequence is provided by the user, set by default as F27.")

  FWD <- "AGRGTTYGATYMTGGCTCAG" # F27
} else {
  print("Forward primer sequence is provided by the user.")

  FWD <- config$primers$fwd
}

print(paste0("Forward primer sequence is ", FWD, "."))

# Reverse primer

if (is.null(config$primers$rev)) {
  print("No reverse primer sequence is provided by the user, set by default as R1492.")

  REV <- "RGYTACCTTGTTACGACTT" # R1492
} else {
  print("Reverse primer sequence is provided by the user.")

  REV <- config$primers$rev
}

print(paste0("Reverse primer sequence is ", REV, "."))

# Setting files in and out paths

folder_deprimed <- file.path(
  output_intermediary,
  "deprimed"
)

deprimed_files <- file.path(
  folder_deprimed,
  "unfiltered",
  paste0(
    samples,
    ".fastq.gz"
  )
)

# Removing the primers

print("Starting depriming")

track <- removePrimers(
  fn = input_files,
  fout = deprimed_files,
  primer.fwd = FWD,
  primer.rev = dada2:::rc(REV), # Being reverse complemented
  compress = TRUE, # Default TRUE
  orient = TRUE, # Default TRUE
  verbose = TRUE # Default False
)

print("Depriming done!")

rownames(track) <- samples

sanity_check(
  files = deprimed_files,
  step_name = "deprimed",
  track = track,
  output_folder = output_intermediary
)

#### --Filtering####

# Setting files in and out paths

filtered_files <- file.path(
  folder_deprimed,
  "filtered",
  paste0(
    samples,
    ".fastq.gz"
  )
)

# Filtering

print("Starting filtering")

track <- filterAndTrim(
  fwd = deprimed_files,
  filt = filtered_files,
  minQ = 3,
  minLen = 1000,
  maxLen = 1600,
  maxN = 0,
  rm.phix = FALSE, # True by default
  maxEE = 2,
  verbose = TRUE,
  multithread = args$threads
)

print("Filtering done!")

rownames(track) <- str_remove(rownames(track), ".fastq.gz")

sanity_check(
  files = filtered_files,
  step_name = "filtered",
  track = track,
  output_folder = output_intermediary
)

#### --Dereplicating, learning errors rates and denoising ####

print("Start dereplicating, learning error rates and denoising")

if (!(args$big_data)) {
  print("Using normal worklow")

  #### --Dereplication####

  print("Starting dereplication")

  derep <- derepFastq(
    fls = filtered_files,
    verbose = TRUE
  )

  names(derep) <- str_remove(names(derep), ".fastq.gz")

  print("Dereplication done!")

  # Saving data and track

  saving_data(
    data = derep,
    data_name = "dereplicated",
    output_folder = output_intermediary,
    seq_col = "uniques",
    track = track
  )

  #### --Learning the error rates####

  print("Starting learning error rates")

  err <- learn_errors_step(
    derep = derep,
    output_folder = output_folder,
    output_intermediary = output_intermediary,
    threads = args$threads
  )

  print("Learning error rates done!")

  #### --Denoising####

  print("Starting denoising")

  dd <- dada(
    derep,
    err = err,
    pool = pooling, # pseudo-pooling https://benjjneb.github.io/dada2/pseudo.html see PSEUDO_PREVALENCE and PSEUDO_ABUNDANCE
    BAND_SIZE = 32,
    multithread = args$threads
  )

  print("Denoising done!")

  # Saving data and track

  saving_data(
    data = dd,
    data_name = "denoised",
    output_folder = output_intermediary,
    seq_col = "denoised",
    track = track
  )
} else {
  print("Using big data worklow")

  #### --Learning the error rates####

  print("Starting learning error rates")

  err <- learn_errors_step(
    derep = filtered_files,
    output_folder = output_folder,
    output_intermediary = output_intermediary,
    threads = args$threads
  )

  print("Learning error rates done!")

  #### --Denoising and dereplication per samples####

  names(filtered_files) <- samples

  derep <- list()

  dd <- list()

  for (sample in samples) {
    cat(paste0("Processing: ", sample, "\n"))

    print("Starting dereplication")

    derep[[sample]] <- derepFastq(
      fls = filtered_files[[sample]],
      verbose = TRUE
    )

    print("Dereplication done!")

    print("Starting denoising")

    dd[[sample]] <- dada(
      derep[[sample]],
      err = err,
      pool = pooling, # pseudo-pooling https://benjjneb.github.io/dada2/pseudo.html see PSEUDO_PREVALENCE and PSEUDO_ABUNDANCE
      BAND_SIZE = 32,
      multithread = args$threads
    )

    print("Denoising done!")
  }

  saving_data(
    data = derep,
    data_name = "dereplicated",
    output_folder = output_intermediary,
    seq_col = "uniques"
  )

  saving_data(
    data = dd,
    data_name = "denoised",
    output_folder = output_intermediary,
    seq_col = "denoised"
  )
}

print("Dereplicating, learning error rates and denoising done!")

#### --Intermediary saving####

print("Saving R workspace...")

save.image(
  file = file.path(
    output_intermediary,
    "read2abd.RData"
  )
)

print("Saving R workspace done!")

#### --Creating abundance table####

print("Creating and saving abundance table...")

abd <- makeSequenceTable(dd)

write.table(
  x = t(abd),
  file = file.path(
    output_intermediary,
    "abd_undechimerized.tsv"
  ),
  sep = "\t",
  row.names = T
)

saveRDS(
  object = abd,
  file = file.path(
    output_intermediary,
    "abd_undechimerized.rds"
  )
)

print("Creating and saving abundance table done!")
