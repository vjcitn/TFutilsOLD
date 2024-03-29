% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grabTab.R
\name{grabTab}
\alias{grabTab}
\title{create table of TF targets and related metadata}
\usage{
grabTab(tfstub = "STAT1", gscoll = TFutils::tftColl,
  orgdb = org.Hs.eg.db::org.Hs.eg.db,
  gwrngs = TFutils::gwascat_hg19_chr17)
}
\arguments{
\item{tfstub}{character(1) gene-like symbol for TF; will be grepped in names(gscoll)}

\item{gscoll}{a GSEABase GeneSetCollection}

\item{orgdb}{an instance of OrgDb as defined in AnnotationDbi}

\item{gwrngs}{a GRanges representing EBI gwascat, must have \code{DISEASE/TRAIT}, \code{MAPPED_GENE}}
}
\value{
data.frame instance
}
\description{
create table of TF targets and related metadata
}
\note{
This function will link together information on targets of a
given TF to the GWAS catalog.
}
\examples{
gt = grabTab("VDR", gscoll=TFutils::tftColl,
   orgdb=org.Hs.eg.db::org.Hs.eg.db, gwrngs=TFutils::gwascat_hg19_chr17)
dim(gt)
head(gt)
}
