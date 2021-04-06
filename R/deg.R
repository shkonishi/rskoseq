#' Differential Expression Gene extraction using TCC package
#' @description Differential expression gene using several methods from the pipeline of TCC packages.
#' @usage deg(dat, gp, method, normalized, fdr)
#' @param dat dataframe: RNA-seq count table.  row: samples, column: genes
#' @param gp factor: group of data. E.g. factor(c(1,1,2,2)); factor(c('A','A','B','B'))
#' @param method integer: Choose from the six pipe line numbers below 1:'DEGES/TbT', 2:'DEGES/edgeR', 3:'iDEGES/edgeR',
#' 4: 'DEGES/DESeq', 5: iDEGES/DESeq, 6: iDEGES/DESeq2. The default value is 2
#' @param normalized logical[defalut FALSE]: if the 'dat' has already normalized set to be TRUE.
#' @param fdr numeric[defalut 0.05]: fdr value.
#' @return TCC object and a result of 'TCC::estimateDE', without normalized count data.
#' @import TCC
#' @importFrom dplyr %>%
#' @importFrom plyr .
#' @examples \dontrun{
#' # sample data of rna-seq
#' fpkm_rep3 <- rskodat::fpkm[1:6]
#' fpkm_rep2 <- rskodat::fpkm[c(1,2,4,5)]
#' fpkm_norep <- rskodat::fpkm[c(1,4)]
#'
#' # aruguments
#' idx1 <- rep(c(1,2), each = 3)
#' idx2 <- rep(c(1,2), each = 2)
#' idx3 <- c(1,2)
#'
#' # two-group with replicate
#' # If you choose 4,5,or 6 which using DESeq, you should load library TCC.
#' res1 <- deg(dat = fpkm_rep3, gp = idx1, fdr = 0.05, method = 1)
#' res2 <- deg(dat = fpkm_rep3, gp = idx1, fdr = 0.05, method = 2)
#' res3 <- deg(dat = fpkm_rep3, gp = idx1, fdr = 0.05, method = 3)
#' res4 <- deg(dat = fpkm_rep3, gp = idx1, fdr = 0.05, method = 4)
#' res5 <- deg(dat = fpkm_rep3, gp = idx1, fdr = 0.05, method = 5)
#' res6 <- deg(dat = fpkm_rep3, gp = idx1, fdr = 0.05, method = 6)
#'
#' # get result of estimateDE
#' head(res3$deg)
#'
#' # draw MA-plot
#' TCC::plot.TCC(res3$tcc)
#'
#' # redraw with ggplot2
#' res3$res[c]
#' deg_dat <- res3m$deg %>% filter(comp %in% levels(.$comp)[1:3])
#' num_de <- res3m$num_deg %>% filter(comp %in% levels(.$comp)[1:3])
#'
#' ggplot2::ggplot(deg_dat, ggplot2::aes(x = a.value, y = m.value, colour = fct)) +
#'   ggplot2::geom_point(size = 0.3) +
#'   ggplot2::scale_color_manual(values =
#'     c("non-DEG"="#BEBEBE80", "up"="#FF000080", "down" ="#0000FF80")) +
#'   ggplot2::theme_minimal(base_size = 15) +
#'   ggplot2::labs(x = "A=(log2(G2)+log2(G1))/2", y = "M=log2(G2)-log2(G1)", colour = "") +
#'   ggplot2::theme(legend.position = "top") +
#'   ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha=1, size = 5))) +
#'   ggplot2::geom_text(data = add_de, ggplot2::aes(x = x, y = y, label = value, colour = key),
#'     size = 5, show.legend = F) +
#'   ggplot2::facet_wrap(~comp, ncol = 3)
#' }
#'
#' @export
deg <- function(dat, gp, method, normalized = F, fdr = 0.05) {
  # create tcc object
  tcc <- methods::new("TCC", dat, gp)

  # if this iput data already normalized
  if (normalized == T) {
    # estimateDE
    if (method == 1) { # 1: 'DEGES/TbT'
      tcc <- TCC::estimateDE(tcc, test.method = "bayseq", FDR = fdr)

    } else if (method == 2 | method == 3) { # 2:'DEGES/edgeR' or 3:'iDEGES/edgeR'
      tcc <- TCC::estimateDE(tcc, test.method = "edger", FDR = fdr)

    } else if (method == 4 | method == 5) { # 4:'DEGES/DESeq' or 5:'iDEGES/DESeq'
      tcc <- TCC::estimateDE(tcc, test.method = "deseq", FDR = fdr)

    } else if (method == 6) { # 6:'iDEGES/DESeq2'
      tcc <- TCC::estimateDE(tcc, test.method = "deseq2", FDR = fdr)
    }

  } else {
    # calcNormFactors & estimateDE
    if (method == 1) { # 1: 'DEGES/TbT'
      samplesize <- 10000
      tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "bayseq",
                                  iteration = 1, samplesize = samplesize)
      tcc <- TCC::estimateDE(tcc, test.method = "bayseq", FDR = fdr)

    } else if (method == 2) { # 2:'DEGES/edgeR'
      tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                                  iteration = 1, FDR = 0.1, floorPDEG = 0.05)
      tcc <- TCC::estimateDE(tcc, test.method = "edger", FDR = fdr)

    } else if (method == 3) { # 3:'iDEGES/edgeR'
      tcc <- TCC::calcNormFactors(tcc, norm.method = "tmm", test.method = "edger",
                                  iteration = 3, FDR = 0.1, floorPDEG = 0.05)
      tcc <- TCC::estimateDE(tcc, test.method = "edger", FDR = fdr)

    } else if (method == 4) { # 4:'DEGES/DESeq'
      tcc <- TCC::calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                                  iteration = 1, FDR = 0.1, floorPDEG = 0.05)
      tcc <- TCC::estimateDE(tcc, test.method = "deseq", FDR = fdr)

    } else if (method == 5) { # 5:'iDEGES/DESeq'
      tcc <- TCC::calcNormFactors(tcc, norm.method = "deseq", test.method = "deseq",
                                  iteration = 3, FDR = 0.1, floorPDEG = 0.05)
      tcc <- TCC::estimateDE(tcc, test.method = "deseq", FDR = fdr)

    } else if (method == 6) { # 5:'iDEGES/DESeq2'
      tcc <- TCC::calcNormFactors(tcc, norm.method = "deseq2", test.method = "deseq2",
                                  iteration = 3, FDR = 0.1, floorPDEG = 0.05)
      tcc <- TCC::estimateDE(tcc, test.method = "deseq2", FDR = fdr)
    }

  }

  # get result of estimateDE
  gene_id <- NULL
  res_de <- TCC::getResult(tcc, sort = TRUE) %>%
    dplyr::mutate(gene_id = as.character(gene_id))
  res_ncnt <- tcc$count
  res <- merge(res_ncnt, res_de, by.x = "row.names", by.y = "gene_id") %>%
    dplyr::arrange(rank)

  return(list(res = res, tcc = tcc))
}
