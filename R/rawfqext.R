#' rawfqext
#' @description Merging fastq files of NEXTSEQ output, lane splitted, created by BaseSpace
#' @usage rawfqext(rawfq_dir, paired, prefix)
#' @param rawfq_dir fastq containing directory
#' @param paired logical, if TRUE is paired end read.
#' @param prefix Characters of sample name, The default value is NULL, collect from directories name.
#' @examples
#' \dontrun{
#' rawfq_dir <- "~/pub/sampledata/fastq/NEXTSEQ"
#' file.remove(list.files(rawfq_dir, "*.gz", full.names = T))
#' prefix <- c("3AD","3AE")
#' rawfqext(rawfq_dir, prefix)
#' }
#' @export

rawfqext <- function(rawfq_dir, paired = FALSE, prefix = NULL){

  # collect fastq containing directory -----
  dirs <- list.files(rawfq_dir,  full.names = T)

  # argument check: prefix----
  if (is.null(prefix)) {
    prefix <- sapply(strsplit(list.files(rawfq_dir), "\\-"), "[", 1)

  } else if (length(dirs) != length(prefix)) {
    stop("length of fasq directory and prefix must to be same")

  }

  # merge files of read per lane -----
  for (i in seq_along(dirs)) {
    if (paired == T) {
      gzr1 <- list.files(dirs[i], ".*_R1_.*.fastq.gz$", full.names = T)
      gzr2 <- list.files(dirs[i], ".*_R2_.*.fastq.gz$", full.names = T)
      r1 <- list.files(dirs[i], ".*_R1_.*.fastq$", full.names = T)
      r2 <- list.files(dirs[i], ".*_R2_.*.fastq$", full.names = T)

      if (!identical(gzr1, character(0)) & identical(r1, character(0))) {
        # command of merge R1 & R2 compressed files ----
        com1g <- paste("cat", paste(gzr1, collapse = " "), ">",
                       paste0(rawfq_dir, "/", prefix[i], "_R1.fastq.gz")
        )
        com2g <- paste("cat", paste(gzr2, collapse = " "), ">",
                       paste0(rawfq_dir, "/", prefix[i], "_R2.fastq.gz")
        )
        # execute command -----
        system(com1g)
        cat(paste0(com1g, "\n"))
        system(com2g)
        cat(paste0(com2g, "\n"))

      } else if (identical(gzr1, character(0)) & !identical(r1, character(0))) {
        # uncompressed R1 & R2 files are merged and compressd to gz file ----
        com1 <- paste("cat", paste(r1, collapse = " "),
                      "| gzip -c >",
                      paste0(rawfq_dir, "/", prefix[i], "_R1.fastq.gz")
        )
        com2 <- paste("cat", paste(r2, collapse = " "),
                      "| gzip -c >",
                      paste0(rawfq_dir, "/", prefix[i], "_R2.fastq.gz")
        )
        system(com1)
        system(com2)

      }

    } else if (paired == F) {
      gzr1 <- list.files(dirs[i], ".*_R1_.*.fastq.gz$", full.names = T)
      r1 <- list.files(dirs[i], ".*_R1_.*.fastq$", full.names = T)

      if (!identical(gzr1, character(0)) & identical(r1, character(0))) {
        # command of merge R1 & R2 compressed files ----
        com1g <- paste("cat", paste(gzr1, collapse = " "), ">",
                       paste0(rawfq_dir, "/", prefix[i], "_R1.fastq.gz")
        )
        # execute command -----
        system(com1g)
        cat(paste0(com1g, "\n"))

      } else if (identical(gzr1, character(0)) & !identical(r1, character(0))) {
        # uncompressed R1 & R2 files are merged and compressd to gz file ----
        com1 <- paste("cat", paste(r1, collapse = " "),
                      "| gzip -c >",
                      paste0(rawfq_dir, "/", prefix[i], "_R1.fastq.gz")
        )
        system(com1)
        cat(paste0(com1, "\n"))

      }

    }

  }
}


