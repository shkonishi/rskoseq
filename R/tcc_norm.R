#' Normarization of RNA-seq count data using tmm
#' @description This functions returns the normarized count data of multiple groups
#' @usage tcc_norm(dat, column, gp, method)
#' @param dat A data frame, matrix. The 'dat' contains columns which names samples and rows contains genes
#' @param column Count data columns without id column. The default is 1:ncol(dat)
#' @param gp replicate group
#' @param method  a pipe line select from   # 1:'DEGES/TbT', 2:'DEGES/edgeR', 3:'iDEGES/edgeR', 4:'DEGES/DESeq',
#' 5:'iDEGES/DESeq', 6:'iDEGES/DESeq2'.
#' @return ggobject which containing result of 'TCC::getNormalizedData' ,
#' @examples
#' \dontrun{
#' # sample data of rna-seq
#' fpkm_rep3 <- rskodat::fpkm
#' fpkm_rep2 <- rskodat::fpkm[c(1,2,4,5)]
#' fpkm_norep <- rskodat::fpkm[c(1, 4)]
#'
#' # arguments
#' gp1 <- rep(1:4, each=3)
#' gp2 <- rep(1:2, each=2)
#' gp3 <- c(1,2)
#'
#' # normalize
#' # if you choose method 4,5, or 6, require the 'TCC' library load.
#' tbt_nfpkm <- rskoseq::tcc_norm(dat = fpkm_rep2, column = 1:ncol(fpkm_rep2), gp = gp2, method = 1)
#' edgr_nfpkm <- rskoseq::tcc_norm(dat = fpkm_rep2, column = 1:ncol(fpkm_rep2), gp = gp2, method = 2)
#' iedgr_nfpkm <- rskoseq::tcc_norm(dat = fpkm_rep2, column = 1:ncol(fpkm_rep2), gp = gp2, method = 3)
#' dsq_nfpkm <- rskoseq::tcc_norm(dat = fpkm_norep, column = 1:ncol(fpkm_norep), gp = gp3, method = 4)
#' idsq_nfpkm <- rskoseq::tcc_norm(dat = fpkm_norep, column = 1:ncol(fpkm_norep), gp = gp3, method = 5)
#'
#' TCC::plot.TCC(TCC::estimateDE(nfpkm1$tcc, test.method = "edger", FDR = 0.05))
#' TCC::plot.TCC(TCC::estimateDE(nfpkm2$tcc, test.method = "edger", FDR = 0.05))
#' TCC::plot.TCC(TCC::estimateDE(nfpkm3$tcc, test.method = "deseq2", FDR = 0.05))
#' }
#'
#' @import TCC
#' @export
tcc_norm <- function(dat, column = 1:ncol(dat), gp, method =2) {

  # argument check: dat ----
  if (is.data.frame(dat)) {
    if (!all(sapply(dat[column], mode) == "numeric")) {
      stop("'dat[column]' just contains count data.")
    } else {
      d <- dat[column]
    }
  } else {
    stop("dat must to be data.frame")
  }

  # argument check: nrow(dat) and length(gp) ----
  if (ncol(d) != length(gp)) {
    stop("'dat' contains columns which names samples, and rows contains genes")
  }

  # create TCC object ----
  tcc <- methods::new("TCC", d, gp)

  # calcNormFactors ----
  if (method == 1) { # 1: 'DEGES/TbT'
    samplesize <- 10000
    tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "bayseq",
                                iteration = 1, samplesize = samplesize)

  } else if (method == 2) { # 2:'DEGES/edgeR'
    tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                                iteration = 1, FDR = 0.1, floorPDEG = 0.05)

  } else if (method == 3) { # 3:'iDEGES/edgeR'
    tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                                iteration = 3, FDR = 0.1, floorPDEG = 0.05)

  } else if (method == 4) { # 4:'DEGES/DESeq'
    tcc <- TCC::calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                                iteration = 1, FDR = 0.1, floorPDEG = 0.05)

  } else if (method == 5) { # 5:'iDEGES/DESeq'
    tcc <- TCC::calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                                iteration = 3, FDR = 0.1, floorPDEG = 0.05)

  } else if (method == 6) { # 6:'iDEGES/DESeq2'
    tcc <- TCC::calcNormFactors(tcc, norm.method = "deseq2", test.method = "deseq2",
                                iteration = 3, FDR = 0.1, floorPDEG = 0.05)
  } else {
    stop("Select from a number of method 1:'DEGES/TbT', 2:'DEGES/edgeR', 3:'iDEGES/edgeR',
         4:'DEGES/DESeq', 5:'iDEGES/DESeq', or 6:'iDEGES/DESeq2' for multi group normarization. ")
  }


  # return getNormalizedData and tcc object ----
  dat[column] <- TCC::getNormalizedData(tcc)
  res <- list(dat, tcc) %>% stats::setNames(., c("ndat", "tcc"))
  return(res)
}
