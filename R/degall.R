#' Differential Expression analysis using all combination
#' @description split count data with index and Differential expression analysis of all comination groups.
#' @usage degall(dat, idx, normalized, meth_norm, param_fdr)
#' @param dat dataframe: RNA-seq count table.  row: samples, column: genes
#' @param idx factor: group of data. E.g. factor(c(1,1,2,2)); factor(c('A','A','B','B'))
#' @param normalized logical: The 'dat' is a normalized count data or not.
#' @param meth_norm integer: Choose from the six pipe line numbers below 1:'DEGES/TbT', 2:'DEGES/edgeR', 3:'iDEGES/edgeR',
#' 4: 'DEGES/DESeq', 5: iDEGES/DESeq, 6: iDEGES/DESeq2. The default value is 2
#' @param param_fdr numeric: fdr value.
#' @return ggplot object which containing result of 'TCC::estimateDE', without normalized count data.
#' @import TCC
#' @importFrom dplyr %>%
#' @importFrom plyr .
#' @examples \dontrun{
#' # sample data of rna-seq
#' fpkm_rep3 <- rskodat::fpkm
#' fpkm_rep2 <- rskodat::fpkm[c(1,2,4,5)]
#' fpkm_norep <- rskodat::fpkm[c(1,4)]
#'
#' # aruguments
#' gp <- sapply(strsplit(names(fpkm_rep3), "_"), function(x) paste(x[1:2], collapse = "_"))
#' index0 <- factor(rep(1:4, each = 3))
#' index1 <- factor(gp, levels=unique(gp))
#' index2 <-factor(rep(1:2, each=2))
#' index3 <- factor(c(1, 2))
#'
#' # two-group with replicate
#' # If you choose 4,5,or 6 which using DESeq, you should load library TCC.
#' res1 <- rskoseq::degall(dat = fpkm_rep2, idx = index2, meth_norm = 1)
#' res2 <- rskoseq::degall(dat = fpkm_rep2, idx = index2, meth_norm = 2)
#' res3 <- rskoseq::degall(dat = fpkm_rep2, idx = index2, meth_norm = 3)
#' res4 <- rskoseq::degall(dat = fpkm_rep2, idx = index2, meth_norm = 4)
#' res5 <- rskoseq::degall(dat = fpkm_rep2, idx = index2, meth_norm = 5)
#' res6 <- rskoseq::degall(dat = fpkm_rep2, idx = index2, meth_norm = 6)
#'
#' # multi-group with replicate
#' res3m <- rskoseq::degall(dat = fpkm_rep3, idx = index0, meth_norm = 3)
#' res6m <- rskoseq::degall(dat = fpkm_rep3, idx = index1, meth_norm = 6)
#'
#' # two samples without replicate
#' res3s <- rskoseq::degall(dat = fpkm_norep, idx = index3, meth_norm = 3)
#' res4s <- rskoseq::degall(dat = fpkm_norep, idx = index3, meth_norm = 4)
#' print(res3s$maplot)
#'
#' # get result of estimateDE
#' head(res3m$deg)
#'
#' # redraw MA-plot
#' library(dplyr)
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
degall <- function(dat, idx, normalized = F, meth_norm = 2, param_fdr = 0.05) {

  # argument check: dat
  if (!is.data.frame(dat) & !is.matrix(dat)) {
    stop("dat is a dataframe or matrix object")
  } else {
    if (ncol(dat) != length(idx)) {
      stop("dat contains columns as samples, and rows as genes.")
    }
  }

  # argument check: meth_norm
  if (!is.null(meth_norm)) {
    if (!any(c(1, 2, 3, 4, 5, 6) %in% meth_norm)) {
      stop("Select a pipeline number froa as follows, \n
           1:'DEGES/TbT', 2:'DEGES/edgeR', 3: 'iDEGES/edger',
           4:'DEGES/DESeq', 5:'iDEGES/DESeq', 6:'iDEGES/DESeq2'.")
    } else {
      cat(paste(c("DEGES/TbT", "DEGES/edgeR", "iDEGES/edgeR", "DEGES/DESeq",
                  "iDEGES/DESeq", "iDEGES/DESeq2")[meth_norm],
                "was selected.\n"))
    }
  }

  # split dataframe using index
  gplist <- split(seq(ncol(dat)), idx)
  dats <- lapply(seq_along(gplist), function(i) {
    tmp <- data.frame(dat[gplist[[i]]]) %>%
      stats::setNames(., paste(levels(idx)[i], 1:length(gplist[[i]]), sep = "_"))
  }) %>%
    stats::setNames(., levels(idx))


  # Round-robin combination
  idx_comb <- utils::combn(levels(idx), 2)
  dat_comb <- lapply(seq(ncol(idx_comb)), function(i) {
    data.frame(dats[[idx_comb[1, i]]], dats[[idx_comb[2, i]]], check.names = F)
  })

  ## combination names
  comb_name <- paste0("G1(", idx_comb[1, ], ") vs G2(", idx_comb[2, ], ")")
  names(dat_comb) <- comb_name

  # all combinatins group
  gmatch <- function(p, v) {
    l <- list()
    for (i in 1:length(p)) {
      l[[i]] <- which(v %in% p[i])
    }
    return(l)
  }
  gps <- vector("list", length = ncol(idx_comb))
  for (i in 1:ncol(idx_comb)) {
    gpi <- gmatch(idx_comb[, i], as.vector(idx))
    gps[[i]] <- unlist(lapply(seq_along(gpi), function(j) {
      rep(j, length(gpi[[j]]))
    }))
  }

  # normalize and DE function
  nde <- function(d, gp, fdr, method) {
    ## create tcc object
    tcc <- methods::new("TCC", d, gp)

    ## calcNormFactors & estimateDE
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

    ## get result of estimateDE
    res_de <- TCC::getResult(tcc, sort = TRUE) %>%
      dplyr::mutate(gene_id = as.character(gene_id))
    res_ncnt <- tcc$count
    res <- merge(res_ncnt, res_de, by.x = "row.names", by.y = "gene_id") %>%
      dplyr::arrange(rank)

    return(res)

  }

  # DE function
  deg <- function(d, gp, fdr, method) {
    tcc <- methods::new("TCC", d, gp)

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

    # get result of estimateDE
    res_de <- TCC::getResult(tcc, sort = TRUE) %>%
      dplyr::mutate(gene_id = as.character(gene_id))
    res_ncnt <- tcc$count
    res <- merge(res_ncnt, res_de, by.x = "row.names", by.y = "gene_id") %>%
      dplyr::arrange(rank)

    return(res)
  }

  # normalization and DE with mapply
  if (normalized == F) {
    # normalize and DE
    res <- mapply(FUN = nde, d = dat_comb, gp = gps,
                  MoreArgs = list(fdr = param_fdr, method = meth_norm), SIMPLIFY = FALSE)

  } else if (normalized == T) {
    res <- mapply(FUN = deg, d = dat_comb, gp = gps,
                  MoreArgs = list(fdr = param_fdr, method = meth_norm), SIMPLIFY = FALSE)
  }

  # extract and merge result of DE
  Row.names <- NULL; gene_id <- NULL; m.value <- NULL; a.value <- NULL;
  q.value <- NULL; estimatedDEG <- NULL; fct <- NULL

  de_list <- lapply(seq_along(res), function(i) {
    res[[i]] %>%
      dplyr::select(c(Row.names, gene_id, m.value, a.value, q.value, rank, estimatedDEG)) %>%
      dplyr::mutate(fct = factor(ifelse(estimatedDEG == 1 & m.value < 0, "down",
                                        ifelse(estimatedDEG == 1, "up", "non-DEG")), levels = c("non-DEG", "up", "down")),
                    comp = factor(rep(comb_name[i], nrow(.))))
  })

  # rbind all deg table
  de_dat <- do.call(rbind, de_list)

  # list of normalized count
  l_ndat <- lapply(seq_along(res), function(i) res[[i]][, 1:(which(names(res[[i]]) == "a.value") - 1)])

  # dataframe for geom_text using facet text position
  x <- NULL; y <- NULL; key <- NULL; value <- NULL
  pos_up <- round(max(de_dat$a.value) - (max(de_dat$a.value) - min(de_dat$a.value)) / 2) - 2
  pos_dw <- round(max(de_dat$a.value) - (max(de_dat$a.value) - min(de_dat$a.value)) / 2) + 2
  pos_y <- ceiling(max(de_dat$m.value))

  de_dat.add <- tapply(de_dat$fct, de_dat$comp, table)
  add_de <- data.frame(comp = names(de_dat.add),
                       up = sapply(de_dat.add, "[", 2),
                       down = sapply(de_dat.add, "[", 3), row.names = NULL) %>%
    tidyr::gather(key = "key", value = "value", -1) %>%
    dplyr::mutate(x = ifelse(key == "up", pos_up, pos_dw), y = rep(pos_y, nrow(.)))


  # ggplot -
    magg <- ggplot2::ggplot(de_dat, ggplot2::aes(x = a.value, y = m.value, colour = fct)) +
      ggplot2::geom_point(size = 0.3) +
      ggplot2::scale_color_manual(values = c("non-DEG" = "#BEBEBE80", "up" = "#FF000080", "down" = "#0000FF80")) +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::labs(x = "A=(log2(G2)+log2(G1))/2", y = "M=log2(G2)-log2(G1)", colour = "") +
      ggplot2::theme(legend.position = "top") +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 5))) +
      ggplot2::geom_text(data = add_de, ggplot2::aes(x = x, y = y, label = value, colour = key),
                         alpha = 1, size = 5, show.legend = F) +
      ggplot2::facet_wrap(~comp, ncol = 3)

  return(list(maplot = magg, deg = de_dat, num_deg = add_de, ndats = l_ndat))
}
