#' Calcuration of N50 and L50 from genome sequence.
#' @description Calculation of N50 and L50 from genomic sequence in fasta format.
#' If the genomic sequence is pseudo chromosomes, L50 is no needed.
#' @usage nx(in_f, N, genome)
#' @param in_f The input file path with fasta format, or DNAStringSet object from Biostrings package.
#' @param N numeric vector: If calcuration for N90(L90), N50(L50), and N10(L10), N = c(90, 50, 10)
#' @param genome numeric: genome size. The default is NULL. The genome size is calculated from the sum of the widths of the input fasta files.
#' @examples
#' # fas <- "~/db/genome/CHOK1GS_HDv1/CHOK1GS_HDv1.dna.toplevel.fa.gz"
#' # rskoseq::nx(in_f = fas, N = 50)
#' @export
nx <- function(in_f, N, genome=NULL){
  # length 'N' must be 1 ----
  if (length(N) != 1) {
    stop("'N' must be a vector of length 1")
  }
  # read fasta file ----
  if (class(in_f) == "DNAStringSet") {
    seq <- in_f
  } else if (file.exists(in_f)) {
    seq <- Biostrings::readDNAStringSet(in_f)
  } else {
    stop("'in_f' must be 'DNAStringSet' object, or fasta file PATH")
  }

  # genome size ----
  if (is.null(genome)) {
    genome <- sum(as.numeric(Biostrings::width(seq)))
  }

  # sort sequence by scaffold length ----
  sortlen <- sort(as.numeric(Biostrings::width(seq)), decreasing = T)

  # Cumulative sum of scaffold length ----
  cmsumlen <- cumsum(sortlen)

  # Scaffold number when cumulative sum reaches N (%) of whole base sequence ----
  LN <- which.min(cmsumlen <= genome/(100/N))

  # Scaffold length when cumulative sum reaches N (%) of whole base sequence ----
  NN <- sortlen[LN]

  # result ----
  res <- stats::setNames(c(NN, LN), c(paste0("N", N), paste0("L", N)))

  # return ----
  return(res)
}
