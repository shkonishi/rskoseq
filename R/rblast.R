#' Local blast search with outfmt 6
#' @description This functions returns a data frame as BLAST output, if multiple databese gave as an arugment, they are all in one,
#'     which fomat is outfmt 6 and additional column. The blast outformat is enforced as follows, '-outfmt "6 std qlen slen sstrand salltitles"'.
#' @usage rblast(in_f, out_f, program, db, ...)
#' @param in_f character: input file path of multifasta format.
#' @param out_f character: output file path or stdout as "-".
#' @param program character: select a blast search program "blastn", "blastp", "blastx", "tblastn", "tblastx"
#' @param db blast database name or file path of fasta: if you gives fasta file path, 'makeblastdb' was executed.
#'  it is not file path. if multiple database using, corresponding out put files path gave as 'out_f'.
#' @param ... additional parameter as character strings E.g. "-num_threads 4 -task megablast".
#' @return blast output of 'outfmt 6' and several other column returnd as named data frame.
#'  if all query "No hits found", could not create data.frame
#' @examples
#' \dontrun{
#' ## create blast data base
#' # com0 <- c("makeblastdb -in TAIR10_chr_all.fas -out TAIR10 -dbtype nucl -parse_seqids")
#' # system(com0, intern = T)
#'
#' # sample fasta of rsko package
#' in_fna <- system.file("extdata", "AtMlos.fna", package="rskodat")
#' bndb <-  "~/db/cdna/TAIR10_cdna"
#' bnout <- rblast(in_f = in_fna, out_f = "-", program = "blastn", db = bndb, "-num_threads 4")
#'
#' # multiple db
#' db1 <- "~/db/cdna/TAIR10_cdna"
#' db2 <- "~/db/cdna/kegg_plant"
#' dbs <- c(db1,db2)
#' in_fna <- system.file("extdata", "AtMlos.fna", package="rskodat")
#' bnouts <- rblast(in_f = in_fna, out_f = "-", program = "blastn", db = dbs, "-num_threads 4")
#' }
#' @export
rblast <- function(in_f, out_f, program, db,  ...){
  # argument check: program PATH ----
  blpath <- suppressWarnings(system(paste("which", program), intern = T))
  if (identical(blpath, character(0))) {
    stop(paste0("There is not ", program, ", or the PATH of this does not found."))
  }

  # argument check: multiple db or not ----
  if (all(out_f != "-")) {
    if (length(out_f) != length(db)) {
      stop("'out_f' must be a vector which has same length of 'db', or '-'. ")
    }
  } else {
    out_f <- rep("-", length(db))
  }

  # data base name
  dbname <- basename(db)

  # blast search
  bout_list <- lapply(seq_along(db), function(i){
    ## system command execution
    com <- paste(blpath, "-query", in_f, "-db", db[i], "-out", out_f[i],
                 "-outfmt", "\'6 std qlen slen sstrand salltitles\'" , ...,  sep = " ")
    cat(paste0(com, " \n"))
    res <- system(com, intern = TRUE)

    ## convert data frame, or read output files
    if (all(out_f == "-")) {
      bout <- lapply(res, function(x){unlist(strsplit(x, "\t"))}) %>%
      {data.frame(do.call(what = rbind, args = .), row.names = NULL, stringsAsFactors = F)} %>%
        transform(., dbname = rep(dbname[i], nrow(.))) %>% # add database name
        dplyr::mutate_at(., c(3, 4, 7, 8, 9, 10, 13, 14), list(~as.numeric)) # convert numeric

    } else {
      bout <- readr::read_delim(out_f[i], delim = "\t", col_names = F, col_types = c("ccnnnnnnnnnnnncc"))
      transform(bout, dbname = rep(dbname[i], nrow(bout)))
    }
  })

  # merge all blast output
  bout <- do.call(rbind, bout_list)

  # add column names
  h1 <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
          "qstart", "qend", "sstart", "send", "evalue", "bitscore") # default outfmt 6 format specifiers
  h2 <- c("qlen", "slen", "sstrand","description") # additional format specifiers
  names(bout)[1:16] <- c(h1,h2)

  # return data frame
  return(bout)
}
