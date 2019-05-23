#' Simple visualization of local blast search output
#' @description The result of blast output using rsko::rblast is dataframe, or
#'     local blast output dat with outfmt '6'. The data has these columns, at least, "qseqid", "sseqid", "pident", "qstart", "qend", "sstart", "send", "qlen", "slen","sstrand".
#' @usage ggblast(blast_out, drawtype)
#' @param blast_out dataframe: result of rskoseq::rblast output
#' @param drawtype character number: 1: facet_wrap with query id, 2: A panel per query, 3: All in one panel
#' @importFrom dplyr %>%
#' @examples \dontrun{
#' ## sample fasta of Biostrings package
#' fp <- system.file("extdata", "AtMlos.fna", package="rskodat")
#' bndb <- "~/db/cdna/TAIR10_cdna" # set you environment
#' bnout <- rskoseq::rblast(fp, "-", "blastn", bndb, "-num_threads 4")
#' res1 <- ggblast(bnout, 1)
#' res2 <- ggblast(bnout, 2)
#' res3 <- ggblast(bnout, 3)
#' }
#' @export
ggblast <- function(blast_out, drawtype){
  # change mode to numeric
  num_col <- c("pident", "length", "qstart", "qend",
               "sstart", "send", "qlen", "slen")

  blast_out <- blast_out %>% dplyr::mutate_at(., dplyr::vars(num_col), dplyr::funs(as.numeric))

  # at least data
  minclmns <- c("qseqid", "sseqid", "pident", "qstart", "qend", "sstart",
               "send", "qlen", "slen", "sstrand", "description")
  if (sum(names(blast_out) %in% minclmns) != 11) {
    stop(paste0("The blast output data must have these columns, which ", "\n",
                paste(minclmns, collapse = ", "), "\n",
                "blast outfmt option should be set like this ", "\n",
                "-outfmt", " \"6 std qlen slen sstrand salltitles\" ")
         )
  }
  blast_out <- blast_out[minclmns]

  # insert a query row to dataframe
  qinpos <- match(unique(blast_out$qseqid), blast_out$qseqid)
  uqdat <- data.frame(unique(blast_out$qseqid), unique(blast_out$qseqid), 100,
                      1, blast_out$qlen[qinpos], 1, blast_out$qlen[qinpos],
                      blast_out$qlen[qinpos], blast_out$qlen[qinpos], "plus",
                      unique(blast_out$qseqid),
                      stringsAsFactors = F) %>%
    stats::setNames(., names(blast_out))

  uqdat_list <- lapply(1:nrow(uqdat), function(i)uqdat[i, ])

  # splited subject data ----
  spdat_list <- lapply(unique(blast_out$qseqid), function(x)blast_out[blast_out$qseqid %in% x,])

  # insert query data to blast result ----
  res <- vector("list", nrow(uqdat))
  for (i in 1:nrow(uqdat)) {res[[i]] <- rbind(uqdat_list[[i]], spdat_list[[i]])}
  qsdat <- Reduce(rbind, res)
  qsdat <- transform(qsdat, ord = 1:nrow(qsdat))

  qseqid = NULL; sseqid = NULL; pident = NULL; sstrand = NULL;
  ord = NULL; qsid = NULL; qstart = NULL; qend = NULL
  ggdat <- qsdat %>%
    dplyr::group_by(qseqid, sseqid) %>%
    dplyr::mutate(rank = dplyr::row_number(sseqid)) %>%
    dplyr::ungroup(qseqid, sseqid) %>%
    dplyr::mutate_at(dplyr::vars(c("pident", "qstart", "qend", "sstart", "send", "qlen", "slen")), as.numeric) %>%
    dplyr::mutate(pident = as.integer(as.character(cut(pident, breaks = c(0,seq(70,100,5)), 8:2, include.lowest = T)))) %>%
    dplyr::mutate(pident = ifelse(sstrand == "plus" & qseqid == sseqid, 1,
                                ifelse(sstrand == "minus", pident + 7, pident))) %>%
    dplyr::mutate(pident = factor(pident),
                  sseqid = factor(sseqid, levels = unique(qsdat$sseqid)),
                  qsid = paste0(qseqid,":", sseqid)) %>%
    dplyr::mutate(qseqid = factor(qseqid, levels = unique(qseqid)))%>%
    dplyr::arrange(dplyr::desc(ord))%>%
    dplyr::mutate(qsid = factor(qsid, levels = unique(qsid)))

  # legend label, color code corresponding to identity ----
  x <- c(0,seq(70,100,5))
  idrange <- rev(sapply(utils::head(seq_along(x), -1), function(i) paste(x[i], x[i+1], sep = ":")))
  vlab <- c("query", idrange, paste0(idrange, "(-)"))
  vcol <- c("gray", rev(RColorBrewer::brewer.pal(7, "Reds")),
            rev(RColorBrewer::brewer.pal(7, "Blues")))
  ## barplot(rep(1, 15), col = vcol, names.arg = vlab, las=2) # check

  if (drawtype == "1") { # ggplot  facet by query columns
    ## color code and legend label
    col_qdat <- vcol[as.integer(levels(factor(ggdat$pident)))]
    lab_qdat <- vlab[as.integer(levels(factor(ggdat$pident)))]
    ## plot
    ggobj1 <- ggplot2::ggplot(data = ggdat, ggplot2::aes(y = qsid, yend = qsid)) +
      ggplot2::theme_bw() +
      ggplot2::geom_segment(data = ggdat, ggplot2::aes(x = qstart, xend = qend, colour = pident), size = 3) +
      ggplot2::scale_color_manual(values = col_qdat, name = "identity",labels = lab_qdat) +
      ggplot2::facet_wrap(~qseqid, ncol = 3, scales = "free") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7),
                     strip.background = ggplot2::element_blank(),
                     strip.text.x = ggplot2::element_blank())
      # ggplot2::scale_y_discrete(labels=function(y)sub(".*:", "", y))
    ## return
    return(ggobj1)

  } else if (drawtype == "2") { # ggplot per query respectively
    ## split blast output by query ----
    splt_ggdat <- lapply(split(rownames(ggdat),ggdat$qseqid),
                         function(x)ggdat[as.numeric(x),])
    ## plot per query ----
    fggdat <- function(qdat){
      col_qdat <- vcol[as.integer(levels(factor(qdat$pident)))]
      lab_qdat <- vlab[as.integer(levels(factor(qdat$pident)))]

      mqdat <- qdat %>%
        dplyr::mutate(sseqid = factor(sseqid, levels = unique(sseqid)))

      ggplot2::ggplot(data = mqdat, ggplot2::aes(y = sseqid, yend = sseqid)) +
        ggplot2::theme_bw() +
        ggplot2::geom_segment(data = mqdat, ggplot2::aes(x = qstart, xend = qend, colour = pident), size = 3) +
        ggplot2::scale_color_manual(values = col_qdat, name = "identity", labels = lab_qdat) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7)) +
        ## draw query start end position
        ggplot2::geom_text(data = mqdat %>% dplyr::filter(rank %% 2 != 0),
                           ggplot2::aes(x = qstart, label = qstart), size = 2,angle = 90,
                           position = ggplot2::position_nudge(0,0.1)) +

        ggplot2::geom_text(data = mqdat %>% dplyr::filter(rank %% 2 != 0),
                           ggplot2::aes(x = qend, label = qend), size = 2,angle = 90,
                           position = ggplot2::position_nudge(0, 0.1)) +

        ggplot2::geom_text(data = mqdat %>% dplyr::filter(rank %% 2 == 0),
                           ggplot2::aes(x = qstart, label = qstart), size = 2,angle = 90,
                           position = ggplot2::position_nudge(0,-0.1)) +

        ggplot2::geom_text(data = mqdat %>% dplyr::filter(rank %% 2 == 0),
                           ggplot2::aes(x = qend, label = qend), size = 2, angle = 90,
                           position = ggplot2::position_nudge(0, -0.1))
    }
    ## return ggplot obj ----
    ggs <- lapply(splt_ggdat, function(x)fggdat(x))
    ggs <- do.call(gridExtra::grid.arrange, c(ggs, list(ncol = 3)))
    return(ggs)
  } else if (drawtype == "3"){
    col_qdat <- vcol[as.integer(levels(factor(ggdat$pident)))]
    lab_qdat <- vlab[as.integer(levels(factor(ggdat$pident)))]
    ## plot ----
    ggobj1 <- ggplot2::ggplot(data = ggdat, ggplot2::aes(y = qsid, yend = qsid)) +
      ggplot2::theme_bw() +
      ggplot2::geom_segment(data = ggdat, ggplot2::aes(x = qstart, xend = qend, colour = pident), size = 3) +
      ggplot2::scale_color_manual(values = col_qdat, name = "identity",labels = lab_qdat)
      # ggplot2::scale_y_discrete(labels=function(y)sub(".*:", "", y))
    ## return ----
    return(ggobj1)
  }

}



