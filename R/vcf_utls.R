#' The utility for handling VCF
#'
#' The utility for handling VCF
#'
#' @name vcf_utls
#'
#' @usage vcf_flt(x, snv, qual)
#' @usage vcf_info(x)
#' @usage vcf_gt(x, labs)
#' @usage vcf_anno(x, cds)
#'
#' @param x data.frame: VCF format
#' @param snv character: "snp", or "indel" [default: snp]
#' @param qual numeric: options of vcf_flt parameters
#' @param labs character: labels of sample name [default: NULL]
#' @param cds DNAStringSet object
#'
#' @examples \dontrun{
#' # sample VCF data
#' vcf <- rskodat::vcf
#'
#' # base filtering of VCF
#' snp <- vcf_flt(vcf, "snp", qual = 900)
#' indel <- vcf_flt(vcf, "indel", qual = 900)
#'
#' # column of 'INFO'
#' info <- vcf_info(snp); dim(info)
#'
#' # column of 'GT' and 'PL'
#' gt <- vcf_gt(snp); dim(gt)
#'
#' # variant annotation
#' ## TransDecoder output fasta
#' orf <- system.file("extdata/longorf.fna.gz", package = "rskodat")
#' longorf <- Biostrings::readDNAStringSet(orf)
#' names(longorf) <- sapply(strsplit(names(longorf), " "), "[", 1)
#'
#' ## annotation of snps
#' res <- vcf_anno(snp, longorf)
#'
#' }
#'


#' @rdname vcf_utls
vcf_flt <- function(x, snv = "snp", qual = 0){

  QUAL <- NULL; FILTER <- NULL; INFO <- NULL
  # snv/indel, quality filtering, extract PASS records,
  if (snv == "snp") {
    x %>%
      dplyr::filter(!grepl("^INDEL", INFO)) %>%
      dplyr::filter(FILTER == "PASS") %>%
      dplyr::filter(QUAL > qual)
  } else if (snv == "indel") {
    x %>%
      dplyr::filter(grepl("^INDEL", INFO)) %>%
      dplyr::filter(FILTER == "PASS") %>%
      dplyr::filter(QUAL > qual)

  }
}

#' @rdname vcf_utls
vcf_info <- function(x){

  if (length(grep("^INDEL", x$INFO)) != 0) {
    print("This VCF containing INDEL.")
    snp <- x[!grepl("^INDEL", x$INFO),]
  } else {
    snp <- x
  }

  # split vector to data.frame
  splt_dat <- function(x, splt, lab = NULL){
    # split the x by 'splt'
    spltx <- strsplit(x, splt)

    # the number of splitted elements
    n_elm <- unique(sapply(spltx, length))

    if (length(n_elm) == 1) { # the length of splitted elements are same
      # convert to data.frame
      d <- data.frame(do.call(rbind, spltx), stringsAsFactors = F)

      # add labels to the data.frame
      if (is.null(lab)) {
        lab <- paste0("v", 1:n_elm)
        stats::setNames(d, lab)
      } else if (length(lab) == n_elm) {
        stats::setNames(d, lab)
      } else {
        stop("the length of 'lab' and splitted elements must to be the same.")
      }

    } else {# different length of splitted elements
      if (is.null(names(spltx))) {
        spltx <- stats::setNames(spltx, 1:length(spltx))
        if (is.null(lab)) {
          utils::stack(spltx)
        } else if (length(lab) == 2) {
          melt_dat <- utils::stack(spltx)
          stats::setNames(melt_dat, lab)
        } else {
          stop("the length of 'lab' must to be the two.")
        }

      } else {
        if (is.null(lab)) {
          utils::stack(spltx)
        } else if (length(lab) == 2) {
          melt_dat <- utils::stack(spltx)
          stats::setNames(melt_dat, lab)
        } else {
          stop("the length of 'lab' must to be the two.")
        }
      }
    }
  }

  # numeric columns
  num_col <- c("QUAL", "AC1", "AF1", "BQB", "DP",
               "FQ", "HWE", "MQ", "MQB", "MQSB", "RPB",
               "SGB", "VDB", "rf", "rr", "af", "ar")

  ind <- NULL
  info <- splt_dat(snp$INFO, ";") %>%
    dplyr::mutate(id = snp$CHROM[ind]) %>%
    {cbind(splt_dat(.$values, "=", c("k","v")), .)} %>%
    dplyr::select(4,5,1,2) %>%
    tidyr::spread(., "k", "v") %>%
    {cbind(., splt_dat(.$DP4, ",", c("rf","rr","af","ar")))}

  cbind(snp[1:7], info[-1:-2]) %>%
    dplyr::mutate_at(num_col, as.numeric)

}

#' @rdname vcf_utls
vcf_gt <- function(x, labs = NULL){
  gtpl <- x[10:ncol(x)] # extract columns of GT & PL
  nsmp <- length(10:ncol(x)) # number of samples
  if (is.null(labs)) {
    labs <- c(paste0("gt.",names(x)[10:ncol(x)]),
              paste0("pl.",names(x)[10:ncol(x)]))
  }

  # split vector to data.frame
  splt_dat <- function(x, splt, lab = NULL){
    # split the x by 'splt'
    spltx <- strsplit(x, splt)

    # the number of splitted elements
    n_elm <- unique(sapply(spltx, length))

    if (length(n_elm) == 1) { # the length of splitted elements are same
      # convert to data.frame
      d <- data.frame(do.call(rbind, spltx), stringsAsFactors = F)

      # add labels to the data.frame
      if (is.null(lab)) {
        lab <- paste0("v", 1:n_elm)
        stats::setNames(d, lab)
      } else if (length(lab) == n_elm) {
        stats::setNames(d, lab)
      } else {
        stop("the length of 'lab' and splitted elements must to be the same.")
      }

    } else {# different length of splitted elements
      if (is.null(names(spltx))) {
        spltx <- stats::setNames(spltx, 1:length(spltx))
        if (is.null(lab)) {
          utils::stack(spltx)
        } else if (length(lab) == 2) {
          melt_dat <- utils::stack(spltx)
          stats::setNames(melt_dat, lab)
        } else {
          stop("the length of 'lab' must to be the two.")
        }

      } else {
        if (is.null(lab)) {
          utils::stack(spltx)
        } else if (length(lab) == 2) {
          melt_dat <- utils::stack(spltx)
          stats::setNames(melt_dat, lab)
        } else {
          stop("the length of 'lab' must to be the two.")
        }
      }
    }
  }

  # reformat GT & PL
  lapply(gtpl, splt_dat, ":", c("GT","PL")) %>%
    {
      cbind(do.call(cbind, lapply(., "[", 1)), do.call(cbind, lapply(., "[", 2)))
    } %>%
    stats::setNames(., labs) %>%
    dplyr::mutate_at(.vars = 1:nsmp, .funs = list(~sub("/", "", .))) %>%
    dplyr::bind_cols(x[1:5], .)
}

#' @rdname vcf_utls
vcf_anno <- function(x, cds){
  POS <- NULL; REF <- NULL; ALT <- NULL
  # extract and translate CDS corresponding to 'CHOROM' from VCF
  cds <- cds[match(x$CHROM, names(cds))]
  trcds <- Biostrings::translate(cds[match(x$CHROM, names(cds))], no.init.codon = T)

  # find mismatch position between character stringsets
  seqdiff <- function(seq1, seq2){
    seq <- strsplit(c(seq1, seq2), split = '')
    mismatches <- which(seq[[1]] != seq[[2]])
    return(mismatches)
  }

  # addition of  annotations
  pos_alt <- x[c(1,2,4,5)]
  snv_anno <- unlist(sapply(1:nrow(pos_alt), function(i){
    # extract vcf data
    idx <- as.character(unlist(pos_alt[i,]))
    id <- idx[1]; pos <- as.integer(idx[2]); ref <- idx[3]; alt <- idx[4]
    # id; pos; ref; alt;pos_alt[i,]

    # REF cds, ALT cds
    ref_cds <- toString(cds[[i]])
    alt_cds <- ref_cds
    if (substr(ref_cds, pos, pos) == ref) {
      substr(alt_cds, pos, pos) <- alt
    } else {
      cat(paste0('The "REF" nuc is different from "', ref, '".'))
    }

    # REF aa, ALT aa
    ref_pep <- toString(trcds[[i]])
    alt_pep <- toString(Biostrings::translate(Biostrings::DNAString(alt_cds),
                                              no.init.codon = T))
    mismatch_pos <- seqdiff(ref_pep, alt_pep)

    # compare REF aa, ALT aa -> synonymous or non-synonymous
    if (length(mismatch_pos)) {
      ref_aa <- substr(ref_pep, mismatch_pos, mismatch_pos)
      alt_aa <- substr(alt_pep, mismatch_pos, mismatch_pos)
      paste0(mismatch_pos,":", ref_aa, "->", alt_aa)

    } else {
      "synonymous"
    }
  }))

  # substitution data and aa seq for provean
  prove <- sapply(strsplit(snv_anno, ":|->"), function(x){
    if (length(x) != 1) {
      paste(x[c(2,1,3)], collapse = "")
    } else {
      "synonymous"
    }
  })
  vtrcds <- sapply(trcds, function(x) toString(x))
  prove.seq <- ifelse(snv_anno != "synonymous", vtrcds, NA)

  snv_flt <- pos_alt %>%
    dplyr::mutate(snv_anno = snv_anno,
           provean = prove,
           prove.seq) %>%
    dplyr::mutate(var.cds = paste0(POS,":",REF,"->",ALT)) %>%
    dplyr::select(1:4,8,5,6,7)
  return(snv_flt)
}
