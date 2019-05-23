#' Filtering and trimming of fastq files using ShortRead package
#' @description filtering and trimming of fastq files and output filtered and compressed fastq files
#' @usage fqfilter (fqd, sffx, outd, qs, nhead, ntail, droplen)
#' @param fqd character: the fully path of a directory of fastq files, or the path of fastq files.
#' @param sffx character: suffix of fastq files. The default value is ".fastq.gz"
#' @param outd output directory
#' @param qs ASCII character corresponding to q-score.
#' @param nhead integer: Remove leading nucleotides.
#' @param ntail integer: Remove trailing nucleotides.
#' @param droplen integer: drop length
#' @examples \dontrun{
#' # output filterd and compressed fastq files ##
#' fqd <- system.file("extdata/E-MTAB-1147", package = "ShortRead")
#' sffx <- ".fastq.gz"
#' outd <- "~/pub/sampledata/rnaseq/project1/qcfq"
#' unlink(outd, recursive = T)
#' fqfilter(fqd=fqd, outd=outd, droplen = 30)
#'
#' # filtered fastq
#' qcfq <- list.files(outd, "qc.fastq.gz", full.names =T)
#' }
#' @export
fqfilter <- function(fqd, sffx = ".fastq.gz", outd, qs = "4", nhead = NULL, ntail = NULL, droplen = NULL){

  # collect fastq files ----
  if (all(dir.exists(fqd))) {
    in_fqs <- list.files(fqd, sffx, full.names = T)
  } else if (all(file.exists(fqd))){
    in_fqs <- fqd
  }

  # path of output files filtered and compressed fastq  ---
  if (is.null(outd)){
    out_fqs <- sprintf("%s", sub(sffx, ".qc.fastq.gz", in_fqs))

  } else {
    if (!file.exists(outd)) dir.create(outd)
    out_fqs <- paste0(outd, "/", sprintf("%s", sub(sffx, ".qc.fastq.gz", basename(in_fqs))))

  }

  # prefix of fastq files (sample name) ----
  fn <- sub(sffx, "", basename(in_fqs))

  # open input stream
  n_dis <- vector("numeric", length = length(in_fqs))
  w_dis <- vector("numeric", length = length(in_fqs))

  for (i in seq_along(in_fqs)){
    stream <- ShortRead::FastqStreamer(in_fqs[i])

    # Function Exit Code
    #on.exit(close(stream))

    ## repeat filtering/trimming
    nfil <- numeric(); wfil <- numeric(); j <- 0
    repeat {
      # input chunk
      fq <- ShortRead::yield(stream)
      if (length(fq) == 0) break
      j <- j+1
      len1  <- length(fq)

      # nFilter(threshold)
      fq <- fq[ShortRead::nFilter()(fq)]
      len2 <- length(fq)
      nfil[j] <- len1 - len2

      # Trimming of leading nucleotids: Remove leading nucleotides. If head 13bp, 'start = 13 +1' .
      if (!is.null(nhead)) fq <- ShortRead::narrow(x = fq, start = nhead + 1 )

      # Trimming of trailing nucleotids: Remove trailing nucleotides. If tail 3bp, 'end = width(fq)-3'.
      if (!is.null(ntail)) fq <- ShortRead::narrow(fq, end = ShortRead::width(fq) - ntail)

      ## Quality trimming: nucleotides q-score less than "4" (phred score 20).
      if (!is.null(qs)) fq <- ShortRead::trimTailw(fq, k = 2, a = qs, halfwidth = 2)

      ## Filtering of a read width: Discard the lead which length less than 'droplen'.
      if (!is.null(droplen)) fq <- fq[ShortRead::width(fq) >= droplen]

      wfil[j] <- len2 - length(fq)

      ## overwite fastq object to 'outd' out put directory.
      ShortRead::writeFastq(object = fq, file = out_fqs[i], mode = "a")
    }
    close(stream)
    n_dis[i] <- sum(nfil)
    w_dis[i] <- sum(wfil)

  }
  res <- data.frame(sample=fn, N_read=n_dis, trimmed_fq=w_dis, check.names = F)
  utils::write.table(res, paste0(outd, "/res_flt.txt"), quote = F, sep = "\t", row.names = F)
  #return(res)
  #print(paste(fn, "N-containing", sum(nfil), "discard", "\n", sep = " "))
  #print(paste(fn, "width of trimmed_fq < 20", sum(wfil), "discard\n",sep = " "))
}





