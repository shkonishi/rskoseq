#' Summary of genomic sequences from fasta file
#' @description Summary of genomic sequences from fasta file
#' @usage genome_summary(in_f, N, genome)
#' @param in_f The input file path with fasta format, or DNAStringSet object from Biostrings package.
#' @param N numeric vector: If calcuration for N90(L90), N50(L50), and N10(L10), N = c(90, 50, 10)
#'     The N50 length is defined as the shortest sequence length at 50 percent of the genome. L50 is the number of contigs whose summed length is N50.
#' @param genome numeric: genome size. The default is NULL. The genome size is calculated from the sum of the widths of the input fasta files.
#' @examples
#' \dontrun{
#' fas <- "~/db/genome/CHOK1GS_HDv1/CHOK1GS_HDv1.dna.toplevel.fa.gz"
#' res <- genome_summary(in_f = fas, N = c(90,50,10))
#'
#' # summary of genomic sequences
#' res$summary
#'
#' # base frequency
#' res$base_freq
#'
#' # disribution of contig lengs
#' hist(log10(res$dist_cntg))
#'
#' }
#' @export
genome_summary <- function(in_f, N, genome = NULL){
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

  # the number of contigs or scaffolds ----
  ncontig <-  length(seq)

  # total bases ----
  wseq <- Biostrings::width(seq)
  bp <- sum(as.numeric(wseq))
  max_cntg <- max(wseq)
  min_cntg <- min(wseq)


  res <- rskoseq::nx(in_f = seq, genome)
  nxbp <- res$bp[res$NGX %in% N]
  lxn <- res$LGX[res$NGX %in% N]


  # GC content(%), N content (%) ----
  frq_base <- apply(Biostrings::alphabetFrequency(seq), 2, sum)
  at <- sum(as.numeric(frq_base[names(frq_base) %in% c("A","T")]))
  gc <- sum(as.numeric(frq_base[names(frq_base) %in% c("C","G")]))
  gc_rate <- gc/(at + gc)
  n_rate <- frq_base[which(names(frq_base) == "N")]/bp

  value <- as.character(c(bp, ncontig, max_cntg, min_cntg, nxbp, lxn,
                      round(n_rate, digits = 3), round(gc_rate, digits = 3)))
  tag <- c("Total sequence length(bp)", "Number of contigs", "Longest contig(bp)", "Shortest contig(bp)",
           paste0("N", N), paste0("L", N), "N", "GC")

  gsum <- data.frame(tag, value)

  return(list(summary = gsum, base_freq = frq_base, dist_cntg = wseq))
}
