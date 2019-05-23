#' Counting kmer of nucleotide sequence
#' @description kmer counting from fastq or fasta format file and return kmer table
#' @usage kmer_strm(x, k, n, resfmt)
#' @return kmer table or kmer count at position
#' @param x An input fastq file path or DNAStrinset object.
#' @param k integer: Length of kmer
#' @param n integer: chunk size for splitting fastq data. The default value is 1L for no splitting.
#' @param resfmt character: "total" as total k-mer, or "pos" as k-mer contents at read position.
#' @importFrom dplyr %>%
#' @importFrom plyr .
#' @importFrom foreach %dopar%
#' @examples
#' \dontrun{
#' # fastq which reads are same length
#' fq <- list.files(system.file("extdata/E-MTAB-1147", package = "ShortRead"), full.names = TRUE)
#'
#' # k-mer table
#' kmer_tab <- rskoseq::kmer_strm(fq[1], 7)
#'
#' # k-mer table with parallel processing
#' kmer_tab <- rskoseq::kmer_strm(fq[1], 7, 5000)
#'
#' # k-mer contents at read position
#' kmer_content <- rskoseq::kmer_strm(fq[1], 7, 5000, "pos")
#'
#' # A DNAStringset object as input which reads are different length
#' fa <- ShortRead::readFastq(fq[1]) %>%
#'   ShortRead::trimTails(., k = 2, a = "4", halfwidth = 2) %>%
#'   ShortRead::sread()
#' res <- rskoseq::kmer_strm(fa, 7, 5000)
#' }
#' @export
kmer_strm <- function(x, k, n = 1L, resfmt = "total"){
  # x is a file path of fastq or DNAStringSet object
  if (class(x) != "DNAStringSet") {
    fa <- Biostrings::readDNAStringSet(x, format = "fastq")
  } else if (class(x) == "DNAStringSet") {
    fa <- x
  } else {
    stop("'x' must be a fastq file path or a DNAStringSet object.")
  }

  # function of k-mer counting from DNAStringset object
  kmer_cnt <- function(fa, k, resfmt = "total"){
    # k-mer table per start position
    wd <- Biostrings::width(fa)
    kmer <- NULL; start <- NULL
    kmer_tab <- unlist(stringr::str_split(toString(fa), ", ")) %>%
      {
        if (length(unique(wd)) != 1) {
          lapply(., function(x){
            st <- 1:(nchar(x) - k + 1)
            ed <- st + k - 1
            res <- stringr::str_sub(x, st, ed)
            mxst <- max(wd) - k + 1
            c(res, rep(NA, mxst - length(res)))
          })
        } else {
          starts <- 1:(unique(wd) - k + 1)
          ends <- starts + k - 1
          lapply(., function(x) stringr::str_sub(x, starts, ends))
        }
      } %>%
      do.call(rbind, .) %>%
      tibble::as_tibble(.name_repair = "minimal") %>%
      stats::setNames(., sprintf(paste0("s%0",nchar(ncol(.)), "d"), 1:ncol(.))) %>%
      tidyr::gather(., "start", "kmer") %>%
      dplyr::group_by(start) %>%
      dplyr::count(., kmer) %>%
      tidyr::spread("start", "n") %>%
      dplyr::filter(., !grepl("N", kmer)) %>%
      replace(., is.na(.), 0) %>%
      {
        if (resfmt == "pos") { # k-mer content at read position
          return(.)
        } else if (resfmt == "total") { # total k-mer table
          {dplyr::mutate(., count = rowSums(.[-1]))} %>%
            dplyr::select(kmer, count) %>%
            dplyr::filter(kmer != "0") %>%
            dplyr::arrange(dplyr::desc(count))
        }
      }
    return(kmer_tab)
  }


  # function of splitting fasta sequences to chunk.
  chunk <- function(x, n) split(x, ceiling(seq_along(x)/n))
  if (as.integer(n) != 1L) {
    chnk_fq <- chunk(fa, n)

    # parallel processing of fasta sequence chunks
    cores <- parallel::detectCores(logical = FALSE)
    cluster <- parallel::makeCluster(cores, 'FORK')
    doParallel::registerDoParallel(cluster)
    res_chnks <- foreach::foreach(fa = chnk_fq) %dopar% {kmer_cnt(fa, k, resfmt)}
    parallel::stopCluster(cluster)

    # merge results of parallel processing
    kmer <- NULL; count <- NULL
    if (resfmt == "total") {
      res <- dplyr::bind_rows(res_chnks) %>%
        dplyr::group_by(kmer) %>%
        dplyr::mutate(count = sum(count)) %>%
        dplyr::slice(1) %>%
        dplyr::arrange(dplyr::desc(count))

    } else if (resfmt == "pos") {
      res <- dplyr::bind_rows(res_chnks) %>%
        dplyr::group_by(kmer) %>%
        dplyr::mutate_at(-1, dplyr::funs(sum)) %>%
        dplyr::slice(1)
    }

  } else {
    res <- kmer_cnt(fa, k, resfmt)
  }

  return(res)
}


