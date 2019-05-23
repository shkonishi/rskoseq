#' Fetch data sets from biomart data base
#' @description Fetch biomart data set
#' @usage fetch_bmt(bmt, host_name, dset, attr)
#' @param bmt character: "ENSEMBL_MART_ENSEMBL", "ENSEMBL_MART_MOUSE", "ENSEMBL_MART_SNP",
#'  "ENSEMBL_MART_FUNCGEN", "plants_mart", "plants_variations", "fungi_mart", "fungi_variations",
#'  "protists_mart", "protists_variations", "metazoa_mart", "metazoa_variations"
#' @param host_name character: set host name from "www.ensembl.org", "plants.ensembl.org",
#' "bacteria.ensembl.org", "fungi.ensembl.org","protists.ensembl.org", "metazoa.ensembl.org".
#' search host names from  http://ensemblgenomes.org/info/data_access. If you can not connect, select the mirror.
#' "asia.ensembl.org".
#' @param dset character:
#' @param attr character: The default values are "chromosome_name","start_position", "end_position",
#' "external_gene_name", "description", and "strand"
#' @examples
#' \dontrun{
#' ## examples ##
#' # return of lists of mart
#' res1 <- rskoseq::fetch_bmt(host_name = "plants.ensembl.org", bmt=NULL)
#'
#' # return of data sets of selected mart
#' res2 <- fetch_bmt(bmt ="plants_mart", host_name = "plants.ensembl.org" )
#'
#' # return of biomart with data set
#' res3 <- fetch_bmt(bmt ="plants_mart", host_name = "plants.ensembl.org", "athaliana_eg_gene")
#' head(res3$dat)
#'
#' # select attibutes
#' attr_al <- biomaRt::listAttributes(res3$mart)
#' atr_1 <- c("chromosome_name","ensembl_transcript_id", "ensembl_gene_id",
#'            "ensembl_peptide_id", "external_gene_name", "entrezgene", "description")
#' atr_2 <- c("chromosome_name","start_position", "end_position",
#'            "external_gene_name", "description", "strand")
#' res4 <- fetch_bmt("plants_mart", "plants.ensembl.org", "athaliana_eg_gene", attr=atr_1)
#' res5 <- fetch_bmt("plants_mart", "plants.ensembl.org", "athaliana_eg_gene", attr=atr_2)
#' head(res4$dat)
#' head(res5$dat)
#'
#' # If list of datasets was already created
#' mart_dat <- rskodat::mart_dat
#' mart_dat[grep("Triticum", mart_dat$description),]
#' args <- unlist(mart_dat[197, 1:3])
#' res6 <- rskoseq::fetch_bmt(args[1], args[2], args[3], attr = atr_1)
#' head(res6$dat)
#'
#' # getSequence This function only works when used with the ENSEMBL_MART_ENSEMBL BioMart.
#' res6 <- fetch_bmt("ENSEMBL_MART_ENSEMBL", "www.ensembl.org", "hsapiens_gene_ensembl" )
#' enst <- c("ENST00000622253", "ENST00000610645", "ENST00000620820")
#' res7 <- biomaRt::getSequence(id=enst, type="ensembl_transcript_id", seqType="cdna",
#' mart=res6$mart)
#' }
#' @export
fetch_bmt <- function(bmt="ENSEMBL_MART_ENSEMBL", host_name="www.ensembl.org", dset=NULL,
                      attr = c("chromosome_name","start_position", "end_position",
                               "external_gene_name", "description", "strand")){
  if (is.null(bmt) & is.null(dset)) {
    marts <- biomaRt::listMarts(host = host_name)
    return(marts)

  } else if (is.null(dset)) {
    mart <- biomaRt::useMart(biomart = bmt, host = host_name)
    res_dsets <- biomaRt::listDatasets(mart, verbose = FALSE)
    return(res_dsets)

  } else if (!is.null(bmt) & !is.null(host_name) & !is.null(bmt)) {
    mart <- biomaRt::useMart(biomart = bmt, dataset = dset,host = host_name)
    # attributes
    dat <- biomaRt::getBM(attributes = attr, mart = mart)

    return(list(mart = mart, dat = dat))
  }

}



