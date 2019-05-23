#' GC content per read
#' @description This functions returns ggplot objects of GC content per read
#' @usage gc_read(fqdir, outdir, suffix, prefix, facet_col)
#' @param fqdir character: The file path of fastq files, or dir path, containing fastq files
#' @param outdir output directory
#' @param suffix A pattern of fastq file suffix. The default is ".fastq.gz"
#' @param prefix A vector of samples name. The default values are names of fastq files containing in fqdir, which substitute 'suffix' character.
#' @param facet_col facet of sequence content and quality score per sample plot.
#' @return  The list of ggplot objects returns and write.table for quality assesment with 'qa'.
#' @importFrom dplyr %>%
#' @examples
#' \dontrun{
#' # arguments
#' fqdir <- system.file("extdata/E-MTAB-1147", package = "ShortRead")
#' suffix <- ".fastq.gz"
#' outdir <- "~/pub/sampledata/rnaseq/project1/qa"
#'
#' # execution
#' res_qa <- rskoseq::gc_read(fqdir, outdir, suffix, facet_col=2)
#' do.call(gridExtra::grid.arrange, c(res_qa, list(ncol=2)))
#'
#' }
#' @export

gc_read <- function(fqdir,
                    outdir,
                    suffix=".fastq.gz",
                    prefix = sub(suffix, "", list.files(fqdir, suffix)),
                    facet_col){
  # file exists or not(full path) ----
  if (!all(file.exists(fqdir))) {
    stop("I cannot find these all files.")
  }
  if (!file.exists(outdir)) {
    dir.create(outdir)
    print(paste( "'", outdir, "'", "was created."))
  }

  # fastq files ----
  fqs <- list.files(fqdir, suffix, full.names = T)
  if (identical(fqs, character(0)) & all(file.exists(fqdir))) fqs <- fqdir

  # function: gc sampler ----
  gc_sampler <- function(file_path, n=NULL){
    if (is.null(n)) n <- 1e+6
    sampler_fq <- ShortRead::FastqSampler(file_path, n)
    fqs <- ShortRead::yield(sampler_fq)
    close(sampler_fq)
    # gc content per read
    v_gc <- Biostrings::letterFrequency(ShortRead::sread(fqs), "GC", as.prob = TRUE) * 100
    return(v_gc)
  }
  # function: cbind.fillna
  cbind.fillna <- function(l, fill){
    mlen <- max(sapply(l, length))
    d <- data.frame(sapply(l,  function(x){ c(x, rep(NA, mlen - length(x)))}), check.names = F)
    return(d)
  }

  # dataframe of GC content per read
  GC = NULL; panels <- NULL
  gc_dat <- lapply(fqs, function(x) gc_sampler(file_path = x)) %>%
    cbind.fillna(.) %>%
    stats::setNames(., prefix) %>%
    tidyr::gather(., key = "sample", value = "GC") %>%
    dplyr::mutate(panels = as.integer(factor(sample, levels = unique(sample)))) %>%
    dplyr::mutate(panels = ifelse(max(panels) > 48, 2, 1))


  # pdf
  if (length(unique(gc_dat$sample)) > 48) {
    grDevices::pdf(paste0(outdir, "/", "gc.pdf"), width = 10, height = 5)
    invisible(lapply(seq_along(unique(gc_dat$panels)), function(i) {
      gc_dat %>% dplyr::filter(panels == i) %>%
        ggplot2::ggplot(., ggplot2::aes(x = GC)) +
        ggplot2::geom_histogram(bins = 20, alpha = 0.5) +
        ggplot2::theme_bw(base_size = 15) +
        ggplot2::labs(x = "GC(%)") +
        ggplot2::facet_wrap(~sample, ncol = 8, nrow = 6) -> gc_gg

      print(gc_gg)
    }))
    grDevices::dev.off()

  } else {
    grDevices::pdf(paste0(outdir, "/", "gc.pdf"), width = 10, height = 5)
    gc_gg <- ggplot2::ggplot(gc_dat, ggplot2::aes(x = GC)) +
      ggplot2::geom_histogram(bins = 20, alpha = 0.5) +
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::labs(x = "GC(%)") +
      ggplot2::facet_wrap(~sample, ncol = facet_col)
    print(gc_gg)
    grDevices::dev.off()
  }

  # return
  return(gc_gg)
}

