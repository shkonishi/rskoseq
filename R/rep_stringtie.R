#' Consecutive execution of Transcript assembly and quantification for RNA-Seq using stringtie
#' @description Consecutive processiong of stringtie for multiple samples, then FPKM and cov data table are created.
#'    Execution following from 'rskoseq::project_rnsq', 'rskoseq::rep_hisat2'.
#' @usage rep_stringtie(bamdir, suffix_bam, guide_gff, res_dir, ...)
#' @param bamdir character: the name of alignment directory, sorted-bam files searched from under this directory.
#'     If designate your own sorted-bam containing directory, give the path of the directory
#'     and 'res_dir' must be directory name under the 'bamdir'.
#' @param suffix_bam character: The default is ".sort.bam".
#' @param guide_gff The file path of guide gff.
#' @param res_dir output directory path, the default is 'paste0(dirname(bamdir), "/res_stringtie")'
#' @param ... additional options of stringtie. E.g. "-e"
#' @examples \dontrun{
#' # create project directory
#' alnd <- "~/pub/sampledata/rnaseq/project1/test1.h2"
#' idx <- "~/db/index/hisat2_idx/CriGri_1.0.ens.add"
#' unlink(alnd, recursive = T)
#' rskoseq::project_rnsq("~/pub/sampledata/rnaseq/project1", "test1.h2", "hisat2")
#' system(paste("tree", alnd))
#'
#' # run rep_hisat2
#' rskoseq::rep_hisat2(alndir = alnd, idx = idx, paired = F)
#' system(paste("tree", alnd))
#'
#' # run rep_stringtie
#' bamdir <- "~/pub/sampledata/rnaseq/project1/test1.h2/res_hisat2"
#' guide <- "~/db/index/hisat2_idx/CriGri_1.0.ens.add.gff"
#' rskoseq::rep_stringtie(bamdir = bamdir, guide_gff = guide)
#' system(paste("tree", alnd))
#' }
#' @export
rep_stringtie <- function(bamdir,
                          suffix_bam=".sort.bam", guide_gff,
                          res_dir=paste0(dirname(bamdir), "/res_stringtie"), ...){
  # argument check: stringtie program PATH ----
  if (!any(grep("stringtie", unlist(strsplit(Sys.getenv("PATH"), ":"))))) {
    stop("There is not stringtie program, or the PATH does not found.")
  }

  # argument check: project directory and project name ----
  if (!file.exists(bamdir)) {
    stop(paste0("\'", bamdir, "\'", " does not found."))
  }

  # argument check: guide_gff ----
  if (!file.exists(guide_gff)) {
    stop(paste0("There is not ", guide_gff))
  }

  # argument check: path of sorted bam files and get all samples name----
  if (file.exists(bamdir)) {
    bamfls <- list.files(bamdir, suffix_bam, full.names = T)
    if (identical(bamfls, character(0))) {
      stop(paste0("There is not '.srot.bam' files in ", bamdir,
                  ", or the suffix of these bam files is different from '",  suffix_bam, "'."))
    }
  } else {
    stop(paste0("There is not ", bamdir))
  }

  # collect sample names from bam files, and create path of result gff files. ----
  smps <- sub(".sort.bam", "",
              sapply(strsplit(bamfls, "/"), function(x) utils::tail(x, 1)))
  res_gff <- paste0(res_dir, "/", "gff", "/", smps, ".gff")

  # stringtie additional options ----
  if (!missing(...)) {
    add_op <- paste0(" ", ..., " ")
  }else{
    add_op <- ""
  }

  # stringtie execution ----
  ## commando log file create
  datestrings <- gsub(":", ".", gsub(" ", "_", date()))
  path_comlog <- paste0(dirname(bamdir), "/", "stringtie_", datestrings, "_log.txt")
  file.create(path_comlog)

  ## open path_comlog of
  con <- file(path_comlog, "a")
  writeLines(date(), con)

  ## detect cores
  cores <- parallel::detectCores()

  ## execute
  writeLines("# assembled transcripts ", con)

  for (i in seq_along(bamfls)) {
    com <- paste("stringtie -p", cores, bamfls[i], "-G", guide_gff, "-o", res_gff[i], add_op, sep = " ")
    cat(paste0(com, " \n"))
    system(com, wait = TRUE)

    ## write command log
    writeLines(com, con)
  }

  # stringtie --merge execution ----
  ## create gff list files (it must be full path)
  resgff <- list.files(paste0(res_dir, "/gff"), ".gff", full.names = T)
  resgff_dat <- data.frame(resgff)
  out_f <- paste0(res_dir, "/resgff.txt")
  utils::write.table(resgff_dat, out_f, sep = "\t", quote = F, row.names = F, col.names = F)

  ## stringtie --merge
  com2 <- paste("stringtie --merge -G", guide_gff, "-o", paste0(res_dir, "/merged.gff"), out_f)
  cat(paste0(com2, " \n"))
  system(com2, wait = T)

  ## over write command log
  writeLines(paste0("\n# stringtie --merge \n", com2), con)

  # stringtie execution using merged gff
  ## define output files
  mgff <- paste0(res_dir, "/", "merged.gff")
  res_mgff <- paste0(res_dir, "/", "mgff", "/", smps, "_m.gff")
  res_ballgown <- paste0(res_dir, "/", "ballgown", "/", smps)
  res_tab <- paste0(res_dir, "/", "tab", "/", smps, ".tab")
  ## replicate execution ----
  writeLines("\n# stringtie -eB", con)
  for (i in seq_along(bamfls)) {
    ## create ballgown directory
    if (!file.exists(res_ballgown[i])) {
      dir.create(res_ballgown[i])
    }
    ## command of 'stringtie -eb'
    com3 <-
      paste("stringtie -p", cores,
             bamfls[i],
             "-e",
             "-A", res_tab[i],
             "-b", res_ballgown[i],
             "-G", mgff,
             "-o", res_mgff[i])
    cat(paste0(com3, " \n"))
    system(com3, wait = T)

    # # write command log ----
    writeLines(com3, con)
  }
  close(con)

  # modified merged gff ----
  dat <- utils::read.table(mgff, header = F, sep = "\t", stringsAsFactors = F)
  v9 <- strsplit(dat$V9, "; ")
  tid <- sub("transcript:", "", sub("transcript_id ", "", sapply(v9, function(x) x[grep("transcript_id", x)])))
  rgid <- sapply(v9, function(x){
    rg <- gsub("ref_gene_id gene:|;$", "", x[grep("^ref_gene_id", x)])
    ifelse(identical(rg, character(0)), NA, rg)
  } )
  gid <- sub("gene_id ", "", sapply(v9, function(x)x[grep("^gene_id", x)]))
  gnm <- sapply(v9, function(x){
    nm <- sub("gene_name ", "", x[grep("^gene_name", x)])
    ifelse(identical(nm, character(0)), NA, nm)
  })
  mmgff <- dat %>% dplyr::mutate(tid, rgid, gid, gnm)

  # create fpkm and cov table from Gene abundances files with the -A option ----
  ## read gene abundances files
  tabs <- lapply(res_tab, function(x){
    utils::read.table(x, sep = "\t", header = T, stringsAsFactors = F)
  })

  gfpkm_list <- lapply(seq_along(tabs), function(i){
    tabs[[i]][c("Gene.ID", "FPKM")] %>%
      stats::setNames(., c("Gene.ID", smps[i])) %>%
      dplyr::arrange(Gene.ID)
    })
  gcov_list <- lapply(seq_along(tabs), function(i){
    tabs[[i]][c("Gene.ID", "Coverage")] %>%
      stats::setNames(., c("Gene.ID", smps[i])) %>%
      dplyr::arrange(Gene.ID)
  })




  # # merge gene abundances files ----
  Gene.ID <- NULL; Gene.Name <- NULL; Ref.Gene.ID <- NULL;  Reference <- NULL;
  Start <- NULL; Strand <- NULL; End <- NULL

  f <- function(x, y)dplyr::full_join(x, y, by = "Gene.ID")
  gfpkm <- Reduce(f, gfpkm_list) %>%
    dplyr::left_join(tabs[[1]][c("Gene.ID", "Gene.Name", "Reference", "Strand", "Start", "End")],
                     ., by = "Gene.ID") %>%
    dplyr::mutate(Gene.ID = sub("gene:", "", Gene.ID)) %>%
    dplyr::mutate(Ref.Gene.ID = replace(Gene.ID, grepl("MSTRG.", Gene.ID),
                                        mmgff$rgid[match(grep("MSTRG.", Gene.ID, value = T), mmgff$gid)])) %>%
    dplyr::mutate(Ref.Gene.ID = ifelse(is.na(Ref.Gene.ID), Gene.ID, Ref.Gene.ID)) %>%
    dplyr::arrange(Reference, Start, Ref.Gene.ID) %>%
    dplyr::select(Ref.Gene.ID, Gene.ID, Gene.Name, Reference, Strand, Start, End, smps)

  gcov <- Reduce(f, gcov_list) %>%
    dplyr::left_join(tabs[[1]][c("Gene.ID", "Gene.Name", "Reference", "Strand", "Start", "End")], ., by = "Gene.ID") %>%
    dplyr::mutate(Gene.ID = sub("gene:", "", Gene.ID)) %>%
    dplyr::mutate(Ref.Gene.ID = replace(Gene.ID,
                               grepl("MSTRG.", Gene.ID),
                               mmgff$rgid[match(grep("MSTRG.", Gene.ID, value = T), mmgff$gid)])) %>%
    dplyr::mutate(Ref.Gene.ID = ifelse(is.na(Ref.Gene.ID), Gene.ID, Ref.Gene.ID)) %>%
    dplyr::mutate(gene_name = mmgff$gnm[match(Ref.Gene.ID, mmgff$rgid)]) %>%
    dplyr::arrange(Reference, Start, Ref.Gene.ID) %>%
    dplyr::select(Ref.Gene.ID, Gene.ID, Gene.Name, Reference, Strand, Start, End, smps)


  # # output file ----
  prjn <- basename(dirname( (dirname(bamdir)) ))
  alnd <- basename(dirname(bamdir))
  gfpkmout <- paste0(res_dir, "/", prjn, "_", alnd, "_g_FPKM.txt")
  gcovout <- paste0(res_dir, "/", prjn, "_", alnd, "_g_cov.txt")

  utils::write.table(gfpkm, gfpkmout, quote = F, sep = "\t", row.names = F, col.names = T)
  utils::write.table(gcov, gcovout, quote = F, sep = "\t", row.names = F, col.names = T)


  # create fpkm and cov table from transcript abundances files with -b option ----
  # # read transcript abundances files ----
  t_dats <- lapply(res_ballgown, function(x){
    utils::read.table(paste(x, "t_data.ctab", sep = "/"), sep = "\t", header = T, stringsAsFactors = F)
  })
  tcov_list <- lapply(t_dats, function(x) x[c("t_name", "gene_id", "gene_name", "cov")])
  tfpkm_list <- lapply(t_dats, function(x) x[c("t_name", "gene_id", "gene_name", "FPKM")])

  invisible(lapply(seq_along(tcov_list), function(i) names(tcov_list[[i]])[[4]] <<- smps[i]))
  invisible(lapply(seq_along(tfpkm_list), function(i) names(tfpkm_list[[i]])[[4]] <<- smps[i]))

  # # merge ----
  t_name <- NULL; gene_id <- NULL; gene_name <- NULL
  f <- function(x, y)dplyr::full_join(x, y, by = c("t_name", "gene_id", "gene_name"))
  tfpkm <- Reduce(f, tfpkm_list) %>%
    dplyr::mutate(t_name = sub("transcript:", "", t_name), gene_id = sub("gene:", "", gene_id)) %>%
    dplyr::mutate(rgid = mmgff$rgid[match(t_name, mmgff$tid)]) %>%
    dplyr::mutate(rgid = ifelse(is.na(rgid), gene_id, rgid)) %>%
    dplyr::select(t_name, gene_id, rgid, gene_name, smps)

  tcov <- Reduce(f, tcov_list) %>%
    dplyr::mutate(t_name = sub("transcript:", "", t_name), gene_id = sub("gene:", "", gene_id)) %>%
    dplyr::mutate(rgid = mmgff$rgid[match(t_name, mmgff$tid)]) %>%
    dplyr::mutate(rgid = ifelse(is.na(rgid), gene_id, rgid)) %>%
    dplyr::select(t_name, gene_id, rgid, gene_name, smps)

  # # output file ----
  prjn <- basename(dirname( (dirname(bamdir)) ))
  alnd <- basename(dirname(bamdir))
  fpkmout <- paste0(res_dir, "/", prjn, "_", alnd, "_t_FPKM.txt")
  covout <- paste0(res_dir, "/", prjn, "_", alnd, "_t_cov.txt")

  utils::write.table(tfpkm, fpkmout, quote = F, sep = "\t", row.names = F, col.names = T)
  utils::write.table(tcov, covout, quote = F, sep = "\t", row.names = F, col.names = T)


}
