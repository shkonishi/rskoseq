#' Reverse and/or Complement of Nucleotide sequence as character string.
#'
#' Reverse and/or Complement of Nucleotide sequence as character string.
#'
#' @usage revcomp(x, strand)
#' @param x character string of nucleotide sequence
#' @param strand character: 'rev','comp', or 'revcomp'
#' @examples
#' seq <- "GCATNGCR"
#' revcomp(seq, "rev")
#' revcomp(seq, "comp")
#' revcomp(seq, "revcomp")
#' @export
revcomp <- function(x, strand = "revcomp"){
  x <- toupper(x)
  if (strand == "revcomp") {
    v <- unlist(strsplit(x, ""))
    rx <- paste(rev(v), collapse = "")
    chartr("ATGCYRKMBVDH", "TACGRYMKVBHD", rx)

  } else if (strand == "rev") {
    v <- unlist(strsplit(x, ""))
    paste(rev(v), collapse = "")

  } else if (strand == "comp") {
    chartr("ATGCYRKMBVDH", "TACGRYMKVBHD", x)
  }

}
