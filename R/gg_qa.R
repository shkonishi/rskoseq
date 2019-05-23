#' Quality assesment with ShortRead package and return ggplot object.
#' @description This functions returns ggplot objects and create report of qa. The qa report is created at above of the fastq directory.
#' @usage gg_qa(fqdir, suffix, prefix, qareport, facet_col, outdir, ow)
#' @param fqdir A vector of file path of fastq files, or dir path, containing fastq files
#' @param suffix A pattern of fastq file suffix. The default is ".fastq.gz"
#' @param prefix A vector of samples name. The default values are names of fastq files containing in fqdir, which substitute 'suffix' character.
#' @param qareport logical. reporting or not. The default value is FALSE.
#' @param facet_col integer. facet of sequence content and quality score per sample plot.
#' @param outdir output directory. The default values is same at fqdir.
#' @param ow logical. It must be TRUE if overwritten. The default value is FALSE.
#' @return  The list of ggplot objects returns and write.table for quality assesment with 'qa'.
#' @importFrom dplyr %>%
#' @examples \dontrun{
#' # arguments
#' p <- system.file("extdata/E-MTAB-1147", package = "ShortRead")
#' sffx <- ".fastq.gz"
#' prfx <- sub(sffx,"",list.files(p, sffx))
#'
#' # execution
#' res_qa <- rskoseq::gg_qa(fqdir=p, suffix=sffx, prefix=prfx, facet_col=2)
#' do.call(gridExtra::grid.arrange, c(res_qa, list(ncol=2)))
#' }
#' @export
gg_qa <- function(fqdir,
                  suffix = ".fastq.gz",
                  prefix = sub(suffix, "", list.files(fqdir, suffix)),
                  qareport = FALSE,
                  facet_col = ceiling(sqrt(length(prefix))),
                  outdir = paste0(dirname(fqdir),"/qa"),
                  ow = FALSE) {

  # argument check: file exists or not(full path)
  if (!all(file.exists(fqdir))) {
    stop("I cannot find these all files.")
  }

  # argument check: output directory exists or not
  if (is.null(outdir)) {
    qadat <- ShortRead::qa(dirPath = fqdir, pattern = suffix)

  } else if (!is.null(outdir) & !file.exists(outdir)) {
    dir.create(outdir)

  } else if (file.exists(outdir) & !identical(list.files(outdir), character(0)) & ow == F) {
    stop(paste("There is some file at ", outdir, "\n",
               " 'ow' must be TRUE if overwritten. "))
  }

  # execute qa and reporting
  if (all(dir.exists(fqdir))) {
    qadat <- ShortRead::qa(dirPath = fqdir, pattern = suffix)
  } else {
    qadat <- ShortRead::qa(dirPath = fqdir)
  }

  if (qareport == TRUE) ShortRead::report(qadat, dest = paste0(outdir, "/report"))


  # initialize of tbl object columns
  read <- NULL; pos <- NULL; Base <- NULL; Count <- NULL; Cycle <- NULL; Qscore <- NULL;
  panels <- NULL;  lane <- NULL; value <- NULL; Score <- NULL; sc_percent <- NULL;

  # lane name
  if (length(rownames(qadat[["readCounts"]])) != length(prefix) | identical(prefix, character(0))) {
    stop("The number of prefix and samples must to be the same.")
  }


  # read count data ----
  readat <- qadat[["readCounts"]] %>%
    tibble::rownames_to_column("sample") %>%
    dplyr::mutate(sample = sub(suffix, "", sample)) %>%
    dplyr::mutate(pos = min(read)/2)

  # sequence content tbl object ----
  scrate <- qadat[["perCycle"]]$baseCall %>%
    dplyr::mutate(lane = sub(suffix, "", lane)) %>%
    dplyr::arrange(Base) %>%
    dplyr::mutate(Base = factor(Base, levels = c("A","T","G","C","N"))) %>%
    dplyr::mutate(lane = sub(suffix, "", lane)) %>%
    dplyr::group_by(Cycle, lane) %>%
    dplyr::mutate(sc_percent = (Count/sum(Count)) * 100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(panels = as.integer(factor(lane, levels = unique(lane)))) %>%
    dplyr::mutate(panels = ifelse(panels > 48, 2, 1))

  # Quality Scores tbl object ----
  qsdat <- qadat[["perCycle"]]$quality %>%
    dplyr::mutate(lane = sub(suffix, "", lane)) %>%
    dplyr::group_by(Cycle, lane) %>%
    dplyr::summarise(median = stats::median(rep(Score, Count)), mean = mean(rep(Score, Count))) %>%
    dplyr::ungroup() %>%
    tidyr::gather(key = "Qscore", value = "value", "median", "mean") %>%
    dplyr::mutate(panels = as.integer(factor(lane, levels = unique(lane)))) %>%
    dplyr::mutate(panels = ifelse(panels > 48, 2, 1))

  # Quality Scores all sample ----
  qsdat2 <- qsdat %>% dplyr::filter(Qscore == "mean")

  # frequentSequences: tbl object: ----
  nReads = NULL; nOccurrences = NULL; nReads_percent = NULL
  freqSeq_top10 <- qadat[["frequentSequences"]] %>%
    dplyr::mutate(lane = sub(suffix, "", lane)) %>%
    dplyr::group_by(lane)

  # frequentSequences: fasta object: top10 ----
  freqSeq_fna <- stats::setNames(Biostrings::DNAStringSet(freqSeq_top10$sequence),
                                 paste(freqSeq_top10$lane, freqSeq_top10$count, sep = "|"))
  Biostrings::writeXStringSet(x = freqSeq_fna, filepath = paste0(outdir, "/freqSeq.fna"))

  # frequentSequences: tbl object sequenceDistribution
  freqSeq_dat <- qadat[["sequenceDistribution"]] %>%
    dplyr::group_by(lane) %>%
    dplyr::mutate(nReads_percent = nReads/sum(nReads)*100) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(lane = sub(suffix, "", lane)) %>%
    dplyr::mutate(panels = as.integer(factor(lane, levels = unique(lane)))) %>%
    dplyr::mutate(panels = ifelse(panels > 48, 2, 1))

  # write.table ----
  df_list <- list(readat, scrate, qsdat, freqSeq_top10)
  fls <- paste0(outdir, c("/readat.txt", "/scrate.txt", "/qsdat.txt", "/freqSeq.txt"))
  invisible(lapply(seq_along(fls), function(i){
    utils::write.table(df_list[[i]], fls[[i]], sep = "\t", quote = F, col.names = T, row.names = F)
  }))



  # ggplot: read count ----
  grDevices::pdf(paste0(outdir, "/", "read.pdf"), width = 10, height = 5)
    read_gg <-
      ggplot2::ggplot(data = readat, ggplot2::aes(x = sample, y = read), fill = sample) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_minimal(base_size = 15) +
      ggplot2::labs(x = "") +
      ggplot2::geom_text(ggplot2::aes(y = pos, label = read),
                         color = "white", angle = 90, size = 3.5) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    print(read_gg)
  grDevices::dev.off()

  # ggplot: sequence content per sample ----
  if (nrow(readat) > 48) {
    grDevices::pdf(paste0(outdir, "/", "sc.pdf"), width = 10, height = 5)
    sc_gg <- lapply(seq_along(unique(scrate$panels)), function(i) {
      scrate %>% dplyr::filter(panels == i) %>%
        ggplot2::ggplot(., ggplot2::aes(x = Cycle, y = sc_percent, group = Base, colour = Base)) +
        ggplot2::geom_line() +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::facet_wrap(~lane, ncol = facet_col) ->  sc_gg
      print(sc_gg)
      sc_gg
    })
    grDevices::dev.off()
  } else {
    grDevices::pdf(paste0(outdir, "/", "sc_s.pdf"), width = 10, height = 5)
      ggplot2::ggplot(scrate, ggplot2::aes(x = Cycle, y = sc_percent, group = Base, colour = Base)) +
        ggplot2::geom_line() +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::facet_wrap(~lane, ncol = facet_col) ->  sc_gg
      print(sc_gg)
    grDevices::dev.off()
  }

  # ggplot: fraquency of 'N' per sample ----
  if (nrow(readat) > 48) {
    grDevices::pdf(paste0(outdir, "/", "sc_N.pdf"), width = 10, height = 5)
    nrgg_smp <- lapply(seq_along(unique(scrate$panels)), function(i) {
        nrgg_smp <- scrate %>%
          dplyr::filter(panels == i, Base == "N") %>%
          ggplot2::ggplot(., ggplot2::aes(x = Cycle, y = sc_percent)) +
          ggplot2::geom_point() +
          ggplot2::theme_bw(base_size = 15) +
          ggplot2::labs(y = "N Content(%)") +
          ggplot2::facet_wrap(~lane, ncol = facet_col)
        print(nrgg_smp)
      })
    grDevices::dev.off()
  } else {
    if (!sum(scrate$Base == "N") == 0) {
      grDevices::pdf(paste0(outdir, "/", "sc_N.pdf"), width = 10, height = 5)
      nrgg_smp <- scrate %>%
        dplyr::filter(Base == "N") %>%
        ggplot2::ggplot(ggplot2::aes(x = Cycle, y = sc_percent)) +
        ggplot2::geom_point() +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::labs(y = "N Content(%)") +
        ggplot2::facet_wrap(~lane, ncol = facet_col)
      print(nrgg_smp)
      grDevices::dev.off()
    }
  }
  # ggplot: sequence content per nuc ----
  grDevices::pdf(paste0(outdir, "/", "sc_b.pdf"), width = 10, height = 5)
    sc_gg2 <-
      ggplot2::ggplot(scrate, ggplot2::aes(x = Cycle, y = sc_percent, colour = Base, group = lane)) +
      ggplot2::geom_line(alpha = 0.5) +
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::labs(y = "Sequence Content(%)") +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
      ggplot2::facet_wrap(~Base, ncol = 5)
    print(sc_gg2)
  grDevices::dev.off()

  # ggplot: Quality Scores per sample ----
  if (nrow(readat) > 48) {
    grDevices::pdf(paste0(outdir, "/", "qs_s.pdf"), width = 10, height = 5)
    qs_gg <- lapply(seq_along(unique(qsdat$panels)), function(i){
      qs_g <- qsdat %>%
        dplyr::filter(panels == i) %>%
        ggplot2::ggplot(., ggplot2::aes(x = Cycle, y = value, group = Qscore, colour = Qscore)) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
        ggplot2::facet_wrap(~lane, ncol = 8, nrow = 6)
        print(qs_g)
        qs_g
    })
    grDevices::dev.off()

  } else {
    grDevices::pdf(paste0(outdir, "/", "qs_s.pdf"), width = 10, height = 5)
      qs_gg <- qsdat %>%
        ggplot2::ggplot(., ggplot2::aes(x = Cycle, y = value, group = Qscore, colour = Qscore)) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
        ggplot2::facet_wrap(~lane, ncol = facet_col)
      print(qs_gg)
    grDevices::dev.off()
  }

  # ggplot: Quality Scores all sample ----
  grDevices::pdf(paste0(outdir, "/", "qs_al.pdf"), width = 10, height = 5)
    gttle <- paste0("All ", length(unique(qsdat$lane)), " samples")
    qs_gg2 <-
      ggplot2::ggplot(qsdat2, ggplot2::aes(x = Cycle, y = value))+
      ggplot2::geom_point(size = 0.7, alpha = 0.5) +
      ggplot2::theme_light(base_size = 15) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(title = gttle, y = "Mean Quality Score")
    print(qs_gg2)
  grDevices::dev.off()

  # freqSeq ggplot per sample ----
  if (nrow(readat) > 48) {
    grDevices::pdf(paste0(outdir, "/", "freqSeq.pdf"), width = 10, height = 5)
      gg_frqseq <- lapply(seq_along(unique(freqSeq_dat$panels)), function(i){
        gg_frq <- freqSeq_dat %>%
          dplyr::filter(panels == i) %>%
          ggplot2::ggplot(., ggplot2::aes(x = nOccurrences, y = nReads_percent)) +
          ggplot2::geom_bar(stat = "identity", fill = "grey50") +
          ggplot2::facet_wrap(~lane, ncol = 8, nrow = 6, scale = "free") +
          ggplot2::labs(y = "nReads(%)") +
          ggplot2::scale_x_continuous(breaks = unique(freqSeq_dat$nOccurrences)) +
          ggplot2::theme_minimal(base_size = 15)
        print(gg_frq)
        gg_frq
      })
    grDevices::dev.off()

  } else {
    grDevices::pdf(paste0(outdir, "/", "freqSeq.pdf"), width = 10, height = 5)
      gg_frqseq <- ggplot2::ggplot(freqSeq_dat, ggplot2::aes(x = nOccurrences, y = nReads_percent)) +
        ggplot2::scale_color_manual(values = "grey50") +
        ggplot2::geom_bar(stat = "identity", fill = "grey50") +
        ggplot2::labs(y = "nReads(%)") +
        ggplot2::scale_x_continuous(breaks = unique(freqSeq_dat$nOccurrences)) +
        ggplot2::theme_minimal(base_size = 15) +
        ggplot2::facet_wrap(~lane, ncol = facet_col, scale = "free")
      print(gg_frqseq)
    grDevices::dev.off()
  }

  # return list of ggplot obj and write.table  list of data frame and ggplot objects ----
  if (exists("nrgg_smp")) {
    gg_list <- list(read = read_gg,
                    seq_cnt_smp = sc_gg,
                    seq_cnt_nuc = sc_gg2,
                    n_gg = nrgg_smp,
                    qsc_smp = qs_gg,
                    qsc_all = qs_gg2,
                    freqSeq = gg_frqseq)
  } else {
    gg_list <- list(read = read_gg,
                    seq_cnt_smp = sc_gg,
                    seq_cnt_nuc = sc_gg2,
                    qsc_smp = qs_gg,
                    qsc_all = qs_gg2,
                    freqSeq = gg_frqseq)
  }

  # return
  return(gg_list)

}
