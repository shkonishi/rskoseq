#' Parse hisat2 or bowtie2 log file
#' @description Parse hisat2 log file and extract mapping rate of all samples.
#' @usage maprate(fp, lab, pair, algnr)
#' @param fp file path of bowtie log
#' @param lab character: file name as sample name.
#' @param pair logical: paired end read or nod. The default value is TRUE.
#' @param algnr character: select an alignr from 'hisat2' or 'bowtie2'. The default value is "bowtie2".
#' @return data frame of mapping rate and ggplot object
#' @examples
#' h2.paired.log <- system.file("extdata", "hisat2_log.txt", package = "rskoseq")
#' label <- c("A","B","C","D")
#' res1 <- maprate(fp = h2.paired.log, lab = label)
#' res2 <- maprate(fp = h2.paired.log, pair = T, lab = label, algnr = "hisat2")
#' @export

maprate <- function(fp, lab = NA, pair = F, algnr = "bowtie2"){

  # argument check: hista2 log file
  if(!file.exists(fp)){
    stop("This file does not exist.")
  }

  # argument check: labels length
  lines <- readLines(fp)
  nsample <- length(grep("overall alignment rate", lines))

  if (!length(lab) == nsample) {
    stop("number of samples and label's length are different")
  } else if (all(is.na(lab))) {
    lab <- seq(nsample)
  } else {
    lab <- lab[order(nchar(lab), lab)]
  }

  # read log file ----
  lines <- gsub("^[ ]+|[ ]+$", "",  lines)
  rn <- sub(" reads; of these:", "", lines[grep("reads; of these", lines)])

  prn_lines <- grep(" were paired; of these:", lines, value = T)
  prn <- as.numeric(sapply(strsplit(prn_lines, " "), "[", 1))
  prn_rate <- as.numeric(sub(".*\\((.*)\\%\\).*", "\\1", prn_lines))

  # paired end ----
  if (!identical(prn_lines, character(0))){
    # unmap ----
    unm_lines <- grep(") aligned concordantly 0 times", lines, value = T)
    unm <- as.numeric(sapply(strsplit(unm_lines, " "), "[", 1))
    unm_rate <- as.numeric(sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(unm_lines, " "), "[", 2)))

    # unique ----
    uqm_lines <- grep(" aligned concordantly exactly 1 time", lines, value = T)
    uqm <- as.numeric(sapply(strsplit(uqm_lines, " "), "[", 1))
    uqm_rate <- as.numeric(sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(uqm_lines, " "), "[", 2)))

    # multiple ----
    mm_lines <- grep(" aligned concordantly >1 times", lines, value = T)
    mm <- as.numeric(sapply(strsplit(mm_lines, " "), "[", 1))
    mm_rate <- as.numeric(sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(mm_lines, " "), "[", 2)))

    # reads of either mapped pairs  ----
    mates <- grep("mates make up the pairs; of these:", lines, value = T)
    n_mates <- as.numeric(sapply(strsplit(mates, " "), "[", 1))

    unm_mates <-  grep("aligned 0 times$", lines, value = T)
    n_unm_mates <- as.numeric(sapply(strsplit(unm_mates, " "), "[", 1))
    unm_mates_rate <- as.numeric(sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(unm_mates, " "), "[", 2)))

    uq_mates <- grep("aligned exactly 1 time$", lines, value = T)
    n_uq_mates <- as.numeric(sapply(strsplit(uq_mates, " "), "[", 1))
    uq_mates_rate <- as.numeric(sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(uq_mates, " "), "[", 2)))

    mlt_mates <- grep("aligned >1 times$", lines, value = T)
    n_mlt_mates <- as.numeric(sapply(strsplit(mlt_mates, " "), "[", 1))
    mlt_mates_rate <- as.numeric(sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(mlt_mates, " "), "[", 2)))

    # data.frame of paired-reads mapping rate ----
    if (algnr == "hisat2"){
      # mapping rate table of hisat2 ----
      mrate <- data.frame(id = lab,
                          num_reads = rn,
                          num_preads = prn,
                          num_unm_reads = unm,
                          num_uq_reads = uqm,
                          num_mp_reads = mm,
                          unmapped = round(unm_rate, 1),
                          multiple = round(mm_rate, 1),
                          unique = round(uqm_rate, 1),

                          num_mates = n_mates,
                          num_unmapped_mates = n_unm_mates,
                          num_uq_mates = n_uq_mates,
                          num_multiple_mates = n_mlt_mates,

                          unmapped_mates = round(unm_mates_rate, 1),
                          multiple_mates = round(mlt_mates_rate, 1),
                          unique_mates = round(uq_mates_rate, 1),
                          stringsAsFactors = F
      )
      # plot ----
      id <- NULL; unmapped <- NULL; multiple <- NULL;  key <- NULL; value <- NULL; pos <- NULL
      unmapped_mates <- NULL; unique_mates <- NULL; multiple_mates <- NULL

      ggmrate <- mrate %>%
        dplyr::select(id, unmapped, multiple, unique) %>%
        tidyr::gather(key="key", value="value", -1) %>%
        dplyr::mutate(key = factor(key, levels=c("unmapped", "multiple", "unique"))) %>%
        dplyr::group_by(id) %>%
        dplyr::mutate(pos = rev(cumsum(rev(value)) - rev(value) + rev(value)/2)) %>%
        ggplot2::ggplot(ggplot2::aes (x = id, y = value, fill = key)) +
        ggplot2::theme_minimal() +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
        #ggplot2::geom_text(ggplot2::aes(label=value), col = "white",size = 3, hjust = 0.5, vjust = 3, position = "stack") +
        ggplot2::geom_text(ggplot2::aes(label=value, y=pos), col = "white") +
        ggplot2::theme(axis.text.x =ggplot2::element_text(angle=90, hjust=1)) +
        ggplot2::labs(x="", y="mapping rate(%)", fill="" )

      ggmate <- mrate %>%
        dplyr::select(id, unmapped_mates, multiple_mates, unique_mates) %>%
        tidyr::gather(key="key", value="value", -1) %>%
        dplyr::mutate(key = factor(key, levels=c("unmapped_mates", "multiple_mates", "unique_mates"))) %>%
        dplyr::group_by(id) %>%
        dplyr::mutate(pos = rev(cumsum(rev(value)) - rev(value) + rev(value)/2)) %>%
        ggplot2::ggplot(ggplot2::aes (x = id, y = value, fill = key)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
        ggplot2::theme_minimal() +
        #ggplot2::geom_text(ggplot2::aes(label=value), col = "white",size = 3, hjust = 0.5, vjust = 3, position = "stack") +
        ggplot2::geom_text(ggplot2::aes(label=value, y=pos), col = "white") +
        ggplot2::theme(axis.text.x =ggplot2::element_text(angle=90, hjust=1)) +
        ggplot2::labs(x="", y="mapping rate(%)", fill="" )

      # return ----
      return(list(mrate, ggmrate, ggmate))



    } else if (algnr == "bowtie2"){
      # mapping rate table of bowtie2 ----
      mrate <- data.frame(id = lab,
                          num_reads = rn,
                          num_preads = prn,
                          num_unm_reads = unm,
                          num_uq_reads = uqm,
                          num_mp_reads = mm,
                          unmapped = round(unm_rate, 1),
                          multiple = round(mm_rate, 1),
                          unique = round(uqm_rate, 1),
                          stringsAsFactors = F)
      # plot ----
      id <- NULL; unmapped <- NULL; multiple <- NULL;  key <- NULL; value <- NULL; pos <- NULL
      unmapped_mates <- NULL; unique_mates <- NULL; multiple_mates <- NULL

      ggmrate <- mrate %>%
        dplyr::select(id, unmapped, multiple, unique) %>%
        tidyr::gather(key="key", value="value", -1) %>%
        dplyr::mutate(key = factor(key, levels=c("unmapped", "multiple", "unique"))) %>%
        dplyr::group_by(id) %>%
        dplyr::mutate(pos = rev(cumsum(rev(value)) - rev(value) + rev(value)/2)) %>%
        ggplot2::ggplot(ggplot2::aes (x = id, y = value, fill = key)) +
        ggplot2::theme_minimal() +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
        #ggplot2::geom_text(ggplot2::aes(label=value), col = "white",size = 3, hjust = 0.5, vjust = 3, position = "stack") +
        ggplot2::geom_text(ggplot2::aes(label=value, y=pos), col = "white") +
        ggplot2::theme(axis.text.x =ggplot2::element_text(angle=90, hjust=1)) +
        ggplot2::labs(x="", y="mapping rate(%)", fill="" )


      # return ----
      return(list(mrate, ggmrate))

    } else {
      stop("You can select alignr from 'hisat2' or 'bowtie2'. ")
    }

  } else { # single
    # unmap ----
    unm_lines <- grep(") aligned 0 times", lines, value = T)
    unm <- sapply(strsplit(unm_lines, " "), "[", 1)
    unm_rate <- sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(unm_lines, " "), "[", 2))

    # unique ----
    uqm_lines <- grep(") aligned exactly 1 time", lines, value = T)
    uqm <- sapply(strsplit(uqm_lines, " "), "[", 1)
    uqm_rate <- sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(uqm_lines, " "), "[", 2))

    # multiple ----
    mm_lines <- grep(") aligned >1 times", lines, value = T)
    mm <- sapply(strsplit(mm_lines, " "), "[", 1)
    mm_rate <- sub("\\((.*)\\%\\)", "\\1", sapply(strsplit(mm_lines, " "), "[", 2))

    # data.frame of single-reads mapping rate ----
    mrate <- data.frame(id=lab,
                        nreads=rn,
                        nunmreads = unm,
                        nuqreads = uqm,
                        nmpreads = mm,
                        unmapped = round(as.numeric(unm_rate), 1),
                        multiple = round(as.numeric(mm_rate), 1),
                        unique = round(as.numeric(uqm_rate), 1),
                        stringsAsFactors = F
    )
    # ----
    id <- NULL; unmapped <- NULL; multiple <- NULL;  key <- NULL; value <- NULL; pos <- NULL
    ggmrate <- mrate %>%
      dplyr::select(id, unmapped, multiple, unique) %>%
      tidyr::gather(key="key", value="value", -1) %>%
      dplyr::mutate(key = factor(key, levels=c("unmapped", "multiple", "unique"))) %>%
      dplyr::group_by(id) %>%
      dplyr::mutate(pos = rev(cumsum(rev(value)) - rev(value) + rev(value)/2)) %>%
      ggplot2::ggplot(ggplot2::aes (x = id, y = value, fill = key)) +
      ggplot2::theme_minimal() +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) +
      #ggplot2::geom_text(ggplot2::aes(label=value), col = "white",size = 3, hjust = 0.5, vjust = 3, position = "stack") +
      ggplot2::geom_text(ggplot2::aes(label=value, y=pos), col = "white") +
      ggplot2::theme(axis.text.x =ggplot2::element_text(angle=90, hjust=1)) +
      ggplot2::labs(x="", y="mapping rate(%)", fill="" )
    return(list(mrate, ggmrate))
  }

}

