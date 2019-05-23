#' Continuous execution of sequence alignment using RSEM
#' @description Description
#' @usage rep_rsem (alndir, fqdir, paired, idx_name, suffix_fq, prefix_fq, ...)
#' @param alndir character: output alignment directory path, created by rskoseq::project_rnsq.
#' @param fqdir character: the fully path of fastq files. The default is 'paste0(dirname(alndir), "/fastq")'.
#'     If this directory containing still analyzed fastq and additional fastq files,
#' @param paired logical: paired or single read. If paired end, the names of fastq file must be "_R1" and "_R2".
#' @param idx_name rsem index file path
#' @param suffix_fq suffix of fastq files. The default value is ".fastq.gz"
#' @param prefix_fq character: The default value is `sub(suffix_fq,"", list.files(fqdir))`
#' @param ... additional rsem-calcurate-expression options. E.g.  "--fragment-length-max"
#' @examples \dontrun{
#' # create index
#' rsem-prepare-reference --bowtie2 ref.fasta refname
#' rsem-prepare-reference --bowtie2 --transcript-to-gene-map list.txt ref.fasta refname
#'
#' # single
#' idx <- "~/db/index/rsem_idx/cge25207.add"
#' alnd <- "~/pub/sampledata/rnaseq/project1/res_single_rsm"
#' unlink(alnd, recursive = T)
#' rskoseq::project_rnsq("~/pub/sampledata/rnaseq/project1", "res_single_rsm", "rsem")
#' rskoseq::rep_rsem(alndir=alnd, idx_name=idx)
#'
#' # single (no project -> you must specifie fastq directory)
#' alnd <- "~/pub/sampledata/rnaseq/res_rsm"
#' fqd <- "~/pub/sampledata/rnaseq/project1/fastq"
#' idx <- "~/db/index/rsem_idx/cge25207.add"
#' unlink(alnd, recursive = T)
#' rskoseq::rep_rsem(alndir=alnd, idx_name=idx, fqdir = fqd, suffix_fq = "_R1.fastq.gz")
#'
#' # paired
#' alnd <- "~/pub/sampledata/rnaseq/project1/res_paired_rsm"
#' fqd <- "~/pub/sampledata/rnaseq/project1/paired_fastq"
#' idx <- "~/db/index/rsem_idx/cge25207.add"
#' unlink(alnd, recursive = T)
#' rskoseq::project_rnsq("~/pub/sampledata/rnaseq/project1", "res_paired_rsm", "rsem")
#' rskoseq::rep_rsem(alndir = alnd, fqdir = fqd, paired = T, idx_name = idx)
#'
#' }
#' @export
rep_rsem <- function(alndir,
                     fqdir = paste(dirname(alndir), "fastq", sep = "/"),
                     paired = FALSE,
                     idx_name,
                     suffix_fq = ".fastq.gz",
                     prefix_fq = sub(suffix_fq, "", list.files(fqdir, suffix_fq)),
                      ...){

  # argument check: collect PATH of fastq files in fastq directory ########################
  # # fastq files or directory exist or not ----
  fqdir <- normalizePath(fqdir)
  path_fq <- list.files(fqdir, suffix_fq, full.names = TRUE)
  if (identical(path_fq, character(0))){
    stop(paste("There is not fastq files in", fqdir, ", or the suffix of fastq is different from", suffix_fq, "."))
  }

  # # Collect path of fastq file ----
  if (paired == T){
    r1fqs <- grep("R1", path_fq, value = T)
    r2fqs <- grep("R2", path_fq, value = T)

  } else {
    r1fqs <- path_fq
  }

  # # prefix of fastq files ----
  if (length(r1fqs) != length(prefix_fq)){
    stop("The length of 'prefix_fq' must be the same as the number of fastq files in 'fqdir'.")
  }

  # # if alignment directory is not exists  ----
  if (!file.exists(alndir)) dir.create(path = alndir, recursive = T)
  alndir <- normalizePath(alndir)
  wd <- getwd() # finally return to this directory
  setwd(alndir) # rsem output directory

  # # log files output under the alignment directory ----
  path_maplog <- paste0(alndir, "/", prefix_fq, "_map_log.txt") # default

  # # index file exists or not ----
  if (!file.exists(paste0(idx_name, ".1.bt2"))){
    stop("Can't find index file.")
  }

  # # rsem-calculate-expression path ----
  rsem <- suppressWarnings(system("which rsem-calculate-expression", intern = T))
  if (!length(rsem)){
    stop("Ther is no rsem-calculate-expression, or PATH environmental variable")
  }

  # # check version ----
  # system(paste0(system("which bowtie2", intern=T), " --version"), intern = T)[[1]]
  # system(paste0(system("which samtools", intern=T), " --version"), intern = T)[[1]]

  # # project name and alignment directory name ----
  alnd <- basename(alndir)
  prjn <- basename(dirname(alndir))

  # # command log file ----
  datestrings <- gsub(":", ".", gsub(" ", "_", date()))
  path_comlog <- list.files(alndir, "log.txt", full.names = T)
  if (identical(path_comlog, character(0))){
    path_comlog <- paste0(alndir, "/", prjn, "_", alnd, "_", datestrings, "_log.txt")
    file.create(path_comlog)
  }

  # # detect cores ----
  cores <- parallel::detectCores()

  # # rsem-calcurate-expression additional options ----
  if (!missing(...)){
    add_op <- paste0(" ", ..., " ")
  }else{
    add_op <- ""
  }

  # rsem execution ############################
  # # open connection of command-log file ----
  con <- file(path_comlog, "a")
  writeLines("# rsem", con)

  # # execute command ----
  if (paired == F){
    for (i in seq_along(r1fqs)){
      com <- paste(rsem, "--bowtie2 --sort-bam-by-coordinate -p",
                   cores,
                   add_op,
                   r1fqs[i], idx_name, prefix_fq[i], ">>", path_maplog[i], "2>&1",
                   sep = " ")
      system(com, wait = T)
      cat(paste0(com, " \n"))
      writeLines(com, con)
    }
  } else {
    for (i in seq_along(r1fqs)){
      com <- paste(rsem, "--bowtie2 --sort-bam-by-coordinate -p",
                   cores,
                   add_op,
                   "--paired-end",
                   r1fqs[i], r2fqs[i],
                   idx_name,
                   prefix_fq[i], ">>", path_maplog[i], "2>&1")
      system(com, wait = T)
      cat(paste0(com, " \n"))
      writeLines(com, con)
    }
  }

  # # close connection of command-log file ----
  close(con)

  # file manipulation ############################
  # # merge results files  ----
  transcript_id <- NULL; gene_id <- NULL;
  gen <- list.files(".", ".genes.results")
  iso <- list.files(".", ".isoforms.results")
  ifpkms <- lapply(seq_along(iso), function(i) {
    readr::read_delim(iso[i], delim = "\t", col_names = T) %>%
      dplyr::select("transcript_id", "gene_id", "FPKM") %>%
      stats::setNames(., c("transcript_id", "gene_id", sub(".isoforms.results", "", iso[i]))) %>%
      dplyr::arrange(transcript_id)
  })
  ifpkm <- data.frame(transcript_id = ifpkms[[1]]$transcript_id,
                      gene_id = ifpkms[[1]]$gene_id,
                      do.call(cbind, lapply(ifpkms, function(x)x[, 3])),
                      check.names = F, stringsAsFactors = F)
  # # ecount ----
  iecnts <- lapply(seq_along(iso), function(i) {
    readr::read_delim(iso[i], delim = "\t", col_names = T) %>%
      dplyr::select("transcript_id", "gene_id", "expected_count") %>%
      stats::setNames(., c("transcript_id", "gene_id", sub(".isoforms.results", "", iso[i]))) %>%
      dplyr::arrange(transcript_id)
  })
  iecnt <- data.frame(transcript_id = iecnts[[1]]$transcript_id,
                      gene_id = iecnts[[1]]$gene_id,
                      do.call(cbind, lapply(iecnts, function(x)x[, 3])),
                      check.names = F, stringsAsFactors = F)

  utils::write.table(iecnt, paste0("./", prjn, ".", alnd, ".iecnt.txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  utils::write.table(ifpkm, paste0("./", prjn, ".", alnd, ".ifpkm.txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)

  # # gfpkm ----
  `transcript_id(s)` <- NULL
  gfpkms <- lapply(seq_along(iso), function(i) {
    readr::read_delim(gen[i], delim = "\t", col_names = T) %>%
      dplyr::select("gene_id", `transcript_id(s)`, "FPKM") %>%
      stats::setNames(., c("gene_id", "transcript_ids", sub(".genes.results", "", gen[i]))) %>%
      dplyr::arrange(gene_id)
  })
  gecnts <- lapply(seq_along(iso), function(i) {
    readr::read_delim(gen[i], delim = "\t", col_names = T) %>%
      dplyr::select("gene_id", `transcript_id(s)`, "expected_count") %>%
      stats::setNames(., c("gene_id", "transcript_ids", sub(".genes.results", "", gen[i]))) %>%
      dplyr::arrange(gene_id)
  })

  gfpkm <- data.frame(gene_id = gfpkms[[1]]$gene_id,
                      transcript_ids = gfpkms[[1]]$transcript_ids,
                      do.call(cbind, lapply(gfpkms, function(x)x[, 3])),
                      check.names = F, stringsAsFactors = F)

  gecnt <- data.frame(gene_id = gecnts[[1]]$gene_id,
                      transcript_ids = gfpkms[[1]]$transcript_ids,
                      do.call(cbind, lapply(gecnts, function(x)x[, 3])),
                      check.names = F, stringsAsFactors = F)

  utils::write.table(gecnt, paste0("./", prjn, ".", alnd, ".gecnt.txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)
  utils::write.table(gfpkm, paste0("./", prjn, ".", alnd, ".gfpkm.txt"),
              sep = "\t", quote = F, col.names = T, row.names = F)


  # result files remove to ./resdir ----
  if(!file.exists("./resdir")) dir.create("./resdir")
  fc1 <- sapply(c(gen, iso), function(x) file.rename(x, paste(alndir, "resdir", x, sep = "/")))
  if (all(fc1)) print("all results files removed to 'resdir'. ")

  # all stat directories removed to ./stats ----
  if(!file.exists("./stats")) dir.create("./stats")
  fc2 <- sapply(list.files(alndir, ".stat"), function(x) file.rename(x, paste(alndir, "stats", x, sep = "/")))
  if (all(fc2)) print("all stat directory removed to 'stats'. ")

  # sorted bam files removed to ./sortedbam ----
  if(!file.exists("./sortbam")) dir.create("./sortbam")
  fc3 <- sapply(list.files(alndir, ".sorted.bam"), function(x)file.rename(x, paste(alndir, "sortbam", x, sep = "/")))
  if (all(fc3)){
    fr3 <- file.remove(list.files(alndir, "\\.bam"))
    if (all(fr3)){
      print("all sorted bam files removed to 'sortbam' directory, and other bam files were deleted. ")
    } else {
      print("Could not all sorted bam files removed")
    }
  }

  # mapping rate ----
  system(paste("cat", paste(path_maplog, collapse = " "), "> map_log.txt"))
  mapr <- rskoseq::maprate(fp = "./map_log.txt", pair = paired, lab = prefix_fq, algnr = "bowtie2")
  utils::write.table(mapr[[1]], "./mapping_rate.txt", sep = "\t", quote = F, col.names = T, row.names = F)
  ggplot2::ggsave(filename = "mapping_rate.pdf", plot = mapr[[2]])

  # restore working directory ----
  setwd(wd)

}
