#' Replicate execution of read alignment using hisat2
#' @description The NGS read alignment using hisat2 for multiple samples.
#'     The input and output directory must be created by 'rskoseq::project_rnsq'.
#' @usage rep_hisat2(alndir, idx, project, fqdir, paired, suffix_fq, ...)
#' @param alndir character: the name of alignment directory. results output
#' @param idx character: the fully path of hisat2 index name.
#' @param project logical:default is TRUE. If does not create project directory using rskoseq::project_rnsq, it is FALSE
#' @param fqdir character: the fully path of fastq files. The default is 'paste0(dirname(alndir), "/fastq")'.
#'     If this directory containing still analyzed fastq and additional fastq files,
#' @param paired logical: paired or single read. If paired end, the names of fastq file must be "_R1" and "_R2".
#' @param suffix_fq character: suffix of fastq files. The default is ".fastq.gz"
#' @param ... additional hisat2 options. E.g.  "--no-spliced-alignment --rna-strandness R"
#' @examples \dontrun{
#' # project TRUE (all result created at directories predetermined by 'project_rnsq')
#' alnd <- "~/pub/sampledata/rnaseq/project1/test1.h2"
#' idx <- "~/db/index/hisat2_idx/CriGri_1.0.ens.add"
#' unlink(alnd, recursive = T)
#' rskoseq::project_rnsq("~/pub/sampledata/rnaseq/project1", "test1.h2", "hisat2")
#' rskoseq::rep_hisat2(alndir = alnd, idx = idx, paired = F)
#' system(paste("tree", alnd))
#'
#' # paired
#' alnd <- "~/pub/sampledata/rnaseq/project1/test2.h2"
#' fqd <- "~/pub/sampledata/rnaseq/project1/paired_fastq"
#' unlink(alnd, recursive = T)
#' rskoseq::project_rnsq("~/pub/sampledata/rnaseq/project1", "test2.h2", "hisat2")
#'
#' adopt <- "--rna-strandness RF --no-softclip --dta"
#' rskoseq::rep_hisat2(alndir = alnd, idx = idx, fqdir = fqd, suffix_fq = "_sub.fastq.gz",
#' paired = T, ... = adopt)
#'
#' # project FALSE (all result create under the alignment directory)
#' alnd <- "~/pub/sampledata/rnaseq/project1"
#' fqd <- "~/pub/sampledata/rnaseq/project1/fastq"
#' rep_hisat2(alndir = alnd, project = F, idx = idx, fqdir = fqd, paired = FALSE)
#' }
#' @export
rep_hisat2 <- function(alndir, idx, project = TRUE,
                       fqdir = paste0(dirname(alndir), "/fastq"),
                       paired = FALSE, suffix_fq=".fastq.gz", ...){

  # system command check: hisat2 and samtools program PATH
  hs2c <- suppressWarnings(system("which hisat2", intern = T))
  samc <- suppressWarnings(system("which samtools", intern = T))
  if (hs2c == 1) {
    stop("There is not hisat2 program, or the PATH does not found.")
  }
  if (samc == 1) {
    stop("There is not samtools program, or the PATH does not found.")
  }

  # argument check: hisat2 index exists or not
  if (!file.exists(paste0(idx, ".1.ht2")) & !file.exists(paste0(idx, ".1.ht2l"))) {
    stop(paste0("Thres is not hisat2 index named as ", idx, " ."))
  }

  # argument check: collect PATH of fastq files in fastq directory
  path_fq <- list.files(fqdir, suffix_fq, full.names = TRUE)
  if (identical(path_fq, character(0))) {
    stop(paste("There is not fastq files in", fqdir, ", or the suffix of fastq is different from", suffix_fq, "."))
  }

  # # collect fastq files pqth
  if (paired == T) {
    r1fqs <- grep(paste0("_R1", suffix_fq), path_fq, value = T)
    r2fqs <- grep(paste0("_R2", suffix_fq), path_fq, value = T)
  } else {
    r1fqs <- path_fq
  }

  # # prefix of fastq files
  if (all(grepl(paste0("_R1", suffix_fq), r1fqs))) {
    prefix <- sub(paste0("_R1", suffix_fq), "",
                  sapply(strsplit(r1fqs, "\\/"), function(x) utils::tail(x, 1)))
  } else {
    prefix <- sub(suffix_fq, "", sapply(strsplit(r1fqs, "/"), function(x) utils::tail(x, 1)))
  }

  # log files output under the alignment directory
  # # if alignment directory is not exists
  if (!file.exists(alndir)) dir.create(path = alndir, recursive = T)

  # # create command log file
  datestrings <- gsub(":", ".", gsub(" ", "_", date()))
  aln <- basename(alndir)
  path_prj <- dirname(alndir)
  prjn <- basename(path_prj)
  path_comlog <- paste0(alndir, "/", prjn, "_", aln,"_", datestrings,"_log.txt")
  file.create(path_comlog)

  # # open connection of command log file
  con <- file(path_comlog, "a")
  writeLines(date(), con)

  # # hisat2 log file(created through hisat2 program)
  com_h2log <- paste0(" 2>> ", alndir, "/", "hisat2_log_", datestrings, ".txt ")

  # # detect cores
  cores <- parallel::detectCores()


  # hisat2 additional options
  if (!missing(...)) {
    add_op <- paste0(" ", ..., " ")
  }else{
    add_op <- ""
  }


  # run hisat2
  # # if project truth, all result directory are still exists .
  # # Otherwise project is FALSE, all result files and directory are under the alignment directory

  if (project == TRUE) {
    # # argument check: alignment directory  exists or not ----
    if (!file.exists(alndir)) {
      stop(paste0(" There is not alignment directory '", alndir, "'. \n"))
    }

    # # bam, h2_log.txt, and met files must be in the same directory. ----
    for (i in seq_along(r1fqs)) {

      # # fail align ----
      failaln <- paste0(alndir, "/", "res_hisat2/failalign", "/", prefix[i], ".failalign.fq.gz")

      # # hisat2 meta file ----
      metfile <- paste0(alndir, "/", "res_hisat2", "/",
                        prefix[i], ".hisat.met.txt ")
      # # samtools sam -> bam -> sort ----
      bamdir <- paste0(alndir, "/", "res_hisat2", "/", prefix[i], ".sort.bam")
      bampfx <- paste0(alndir, "/", "res_hisat2", "/", prefix[i], "_sort")
      com_smtools <- paste(samc, "sort -O bam -o", bamdir,"-T", bampfx, "-@",cores )

      # # execute command
      if (paired == TRUE) { # # paired end
        com <- paste0(hs2c, " -p ", cores, " -x ", idx,
                      " --un-conc-gz ", failaln,
                      " --met-file ", metfile,
                      " --dta",
                      add_op,
                      " -1 ", r1fqs[i],
                      " -2 ", r2fqs[i],
                      com_h2log, " | ",
                      com_smtools)
        cat(paste0(com, " \n"))
        return_com <- system(com, intern = T)

        # # if only hisat2 result has statement 1, I don't know whether error occurance. ----
        if (!identical(return_com, character(0))) {
          stop("hisat2 or returned error.  ")
        }

        # # write command to command_log.txt ----
        writeLines(com, con)

      } else if (paired == FALSE) { # single
        com <- paste(hs2c, "-p", cores, "-x", idx,
                      "--un-conc-gz", failaln,
                      "--met-file", metfile,
                      add_op,
                      "-U", r1fqs[i],
                      com_h2log, "|",
                      com_smtools)
        cat(paste0(com, " \n"))
        return_com <- system(com, intern = T)
        # # if only hisat2 result has statement 1, I don't know whether error occurance. ----
        if (!identical(return_com, character(0))) {
          stop("hisat2 or returned error.  ")
        }

        # # write command to command_log.txt ----
        writeLines(com, con)
      }
    }
  } else if (project == FALSE) {
    # # create failalign dir ----
    dir.create(path = paste0(alndir, "/failalign"), recursive = T)

    # # bam, h2_log.txt, and met files must be in the same directory. ----
    for (i in seq_along(r1fqs)) {
      # # fail align ----
      failaln <- paste0(alndir, "/failalign/", prefix[i], ".failalign.fq.gz")

      # # hisat2 meta file ----
      metfile <- paste0(alndir, "/", prefix[i], ".hisat.met.txt ")

      # # samtools sam -> bam -> sort ----
      bamdir <- paste0(alndir, "/", prefix[i], ".sort.bam")
      bampfx <- paste0(alndir, "/", prefix[i], "_sort")
      com_smtools <- paste(samc, "sort -O bam -o", bamdir,"-T", bampfx, "-@",cores )

      # execute command ----
      if (paired == TRUE) { # # paired end
        # # command text ----
        com <- paste0(hs2c, " -p ", cores, " -x ", idx,
                      " --un-conc-gz ", failaln,
                      " --met-file ", metfile,
                      " --dta",
                      add_op,
                      " -1 ", r1fqs[i],
                      " -2 ", r2fqs[i],
                      com_h2log, " | ",
                      com_smtools)
        cat(paste0(com, " \n"))
        return_com <- system(com, intern = T)

        # # if only hisat2 result has statement 1, I don't know whether error occurance. ----
        if (!identical(return_com, character(0))){
          stop("hisat2 or returned error.  ")
        }

        # # write command to command_log.txt ----
        writeLines(com, con)

      } else if (paired == FALSE) { # single
        # # command text ----
        com <- paste0(hs2c, " -p ", cores, " -x ", idx,
                      " --un-conc-gz ", failaln,
                      " --met-file ", metfile,
                      add_op,
                      "-U ", r1fqs[i],
                      com_h2log, " | ",
                      com_smtools)
        cat(paste0(com, " \n"))
        return_com <- system(com, intern = T)
        # # if only hisat2 result has statement 1, I don't know whether error occurance. ----
        if (!identical(return_com, character(0))) {
          stop("hisat2 or returned error.  ")
        }

        # # write command to command_log.txt ----
        writeLines(com, con)
      }
    }

  } else {
    stop("fail")
  }
  close(con)
}
