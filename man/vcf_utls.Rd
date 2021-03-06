% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vcf_utls.R
\name{vcf_utls}
\alias{vcf_utls}
\alias{vcf_flt}
\alias{vcf_info}
\alias{vcf_gt}
\alias{vcf_anno}
\title{The utility for handling VCF}
\usage{
vcf_flt(x, snv, qual)

vcf_info(x)

vcf_gt(x, labs)

vcf_anno(x, cds)

vcf_info(x)

vcf_gt(x, labs = NULL)

vcf_anno(x, cds)
}
\arguments{
\item{x}{data.frame: VCF format}

\item{snv}{character: "snp", or "indel" [default: snp]}

\item{qual}{numeric: options of vcf_flt parameters}

\item{labs}{character: labels of sample name [default: NULL]}

\item{cds}{DNAStringSet object}
}
\description{
The utility for handling VCF
}
\examples{
\dontrun{
# sample VCF data
vcf <- rskodat::vcf

# base filtering of VCF
snp <- vcf_flt(vcf, "snp", qual = 900)
indel <- vcf_flt(vcf, "indel", qual = 900)

# column of 'INFO'
info <- vcf_info(snp); dim(info)

# column of 'GT' and 'PL'
gt <- vcf_gt(snp); dim(gt)

# variant annotation
## TransDecoder output fasta
orf <- system.file("extdata/longorf.fna.gz", package = "rskodat")
longorf <- Biostrings::readDNAStringSet(orf)
names(longorf) <- sapply(strsplit(names(longorf), " "), "[", 1)

## annotation of snps
res <- vcf_anno(snp, longorf)

}

}
