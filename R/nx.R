#' Calcuration of NX and LX from genome sequence.
#' @description Calculation of NX and LX from genomic sequence in fasta format.
#' If the genomic sequence is pseudo chromosomes, L50 is no needed.
#' @usage nx(in_f, genome)
#' @param in_f The input file path with fasta format, or DNAStringSet object from Biostrings package.
#' @param genome numeric: genome size. The default is NULL. The genome size is calculated from the sum of the widths of the input fasta files.
#' @examples \dontrun{
#' fas <- "~/db/genome/CHOK1GS_HDv1/CHOK1GS_HDv1.dna.toplevel.fa.gz"
#' nxdat <- rskoseq::nx(in_f = fas)
#' plot(x = nxdat$NGX, y = nxdat$Mbp, xlab = "NX", ylab = "contig length(Mbp)")
#' abline(v = 50, lty = 3)
#' }
#' @export
nx <- function(in_f, genome = NULL){

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
  sfn <- sapply(1:99, function(n) which.min(cmsumlen <= genome/(100/n)))

  # Scaffold length when cumulative sum reaches N (%) of whole base sequence ----
  sfl <- sortlen[sfn];
  ngdat <- data.frame(NGX = 1:99, LGX = sfn, bp = sfl, Mbp = sfl/10^6)

  # return ----
  return(ngdat)
}
