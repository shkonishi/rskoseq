#' Create directory for RNA-seq pipeline
#' @description Create Project directories for RNA-seq pipeline using hisat2-stringtie.
#'    Then fastq files move to './project/fastq'　directory, followning execution of 'rskoseq::rep_hisat2' or 'rskoseq::rep_rsem', for read mapping.
#' @usage project_rnsq(prjd, alnd, algnr)
#' @param prjd character: path of a project directory name
#' @param alnd character: path of a alignment directory name, the default is "alignment1"
#' @param algnr character: select an alignment program from "hisat2"(default), "rsem", "others"
#' @examples
#' \dontrun{
#' # create new project
#' prj <- "~/pub/sampledata/rnaseq/project1"
#' project_rnsq(prjd = prj, alnd = "test.h2", algnr = "hisat2")
#' project_rnsq(prjd = prj, alnd = "test.rsm", algnr = "rsem")
#'
#' # now there is a project directory, create another alignment directory.
#' system(paste("tree", prj))
#'
#' }
#' @export
project_rnsq <- function(prjd, alnd="alignment1", algnr="hisat2"){
  # argument check: prjd
  prjd <- normalizePath(prjd)
  prjn <- sapply(strsplit(prjd, "/"), function(x) utils::tail(x,1)) # get project name from project path

  # still exist project directory, and create another alignment ----
  if (file.exists(prjd) & is.na(match(alnd, list.files(prjd)))) {
    if (algnr == "hisat2") {
      ## res_hisat2 dir ----
      h2dir <- paste0(prjd, "/", alnd, "/res_hisat2/failalign")
      dir.create(path = h2dir, recursive = T)

      ## res_stringtie dir ----
      bgdir <- paste0(prjd, "/", alnd, "/res_stringtie/ballgown")
      tabdir <- paste0(prjd, "/", alnd, "/res_stringtie/tab")
      gffdir <- paste0(prjd, "/", alnd, "/res_stringtie/gff")
      mgffdir <- paste0(prjd, "/", alnd, "/res_stringtie/mgff")

      dir.create(path = bgdir, recursive = T)
      dir.create(path = tabdir, recursive = T)
      dir.create(path = gffdir, recursive = T)
      dir.create(path = mgffdir, recursive = T)

      # ## command log file ----
      # logfile <- paste0(prjd, "/",alnd, "/", prjn, "_", alnd, "_", "log.txt")
      # file.create(logfile)

    }else if (algnr == "rsem") {
      ## res_rsem, stats, results, bam----
      rsmdir <- paste0(prjd, "/", alnd)
      statsdir <- paste0(prjd, "/", alnd, "/stats")
      resdir <- paste0(prjd, "/", alnd, "/resdir")
      bamdir <- paste0(prjd, "/", alnd, "/sortbam")

      dir.create(path = rsmdir, recursive = T)
      dir.create(path = statsdir, recursive = T)
      dir.create(path = resdir, recursive = T)
      dir.create(path = bamdir, recursive = T)

      # ## command log file ----
      # logfile <- paste0(prjd, "/",alnd, "/", prjn, "_", alnd, "_", "log.txt")
      # file.create(logfile)

    } else if (algnr == "others") {

    }
  # still exist alignment directory, different alignment directory create
  } else if (file.exists(prjd) & !is.na(match(alnd, list.files(prjd)))) {
    stop(paste0("The directory ", prjd, " still exists, and give different name of alignment directory from '", alnd, "'."))


  # create newly project
  } else if (!file.exists(prjd)) {
    ## project directory
    dir.create(prjd)
    dir.create(paste(prjd, "fastq", sep = "/"))
    dir.create(paste(prjd, "qa", sep = "/"))

    if (algnr == "hisat2") {
      ## res_hisat2 dir
      h2dir <- paste0(prjd, "/", alnd, "/res_hisat2/failalign")
      dir.create(path = h2dir, recursive = T)

      ## res_stringtie dir
      bgdir <- paste0(prjd, "/", alnd, "/res_stringtie/ballgown")
      tabdir <- paste0(prjd, "/", alnd, "/res_stringtie/tab")
      gffdir <- paste0(prjd, "/", alnd, "/res_stringtie/gff")
      mgffdir <- paste0(prjd, "/", alnd, "/res_stringtie/mgff")

      dir.create(path = bgdir, recursive = T)
      dir.create(path = tabdir, recursive = T)
      dir.create(path = gffdir, recursive = T)
      dir.create(path = mgffdir, recursive = T)

      # # command log file ----
      # logfile <- paste0(prjd, "/",alnd, "/", prjn, "_", alnd, "_", "log.txt")
      # file.create(logfile)

    } else if (algnr == "rsem"){
      ## res_rsem, stats, results, bam----
      statsdir <- paste0(prjd, "/", alnd, "/stats")
      resdir <- paste0(prjd, "/", alnd, "/resdir")
      bamdir <- paste0(prjd, "/", alnd, "/sortbam")

      dir.create(path = statsdir, recursive = T)
      dir.create(path = resdir, recursive = T)
      dir.create(path = bamdir, recursive = T)

    } else if (algnr == "others"){


    }
  }
}







