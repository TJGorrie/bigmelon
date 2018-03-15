\name{estimateCellCounts}
\alias{estimateCellCounts.gds}
\alias{getquantiles}
\title{
Cell Proportion Estimation using bigmelon
}
\description{
Estimates relative proprotion of pure cell types within a sample, identical to \code{\link[minfi]{estimateCellCounts}}. Currently, only a reference data-set exists for 450k arrays. As a result, if performed on EPIC data, function will convert gds to 450k array dimensions (this will not be memory efficient).
}
\usage{
    estimateCellCounts.gds(
        gds,
        gdPlatform = c("450k", "EPIC", "27k"),
        mn = NULL,
        un = NULL,
        bn = NULL,
        perc = 1,
        compositeCellType = "Blood",
        probeSelect = "auto",
        cellTypes = c("CD8T","CD4T","NK","Bcell","Mono","Gran"),
        referencePlatform = c("IlluminaHumanMethylation450k",
            "IlluminaHumanMethylationEPIC",
            "IlluminaHumanMethylation27k"),
        returnAll = FALSE,
        meanPlot = FALSE,
        verbose=TRUE,
        ...)
}
\arguments{
\item{gds}{An object of class gds.class, which contains (un)normalised methylated and unmethylated intensities}
\item{gdPlatform}{Which micro-array platform was used to analysed samples}
\item{mn}{
'Name' of gdsn node within gds that contains methylated intensities, if NULL it will default to 'methylated' or 'mnsrank' if \code{\link{dasenrank}} was used prior
}
\item{un}{'Name' of gdsn node within gds that contains unmethylated intensities, if NULL it will default to 'unmethylated' or 'unsrank' if \code{\link{dasenrank}} was used prior}

\item{bn}{'Name' of gdsn node within gds that contains un(normalised) beta intensities. If NULL - function will default to 'betas'.}
\item{perc}{Percentage of query-samples to use to normalise reference dataset. This should be 1 unless using a very large data-set which will allow for an increase in performance}
\item{compositeCellType}{Which composite cell type is being deconvoluted. Should be either "Blood", "CordBlood", or "DLPFC"}
  \item{probeSelect}{How should probes be selected to distinguish cell types? Options include
    "both", which selects an equal number (50) of probes (with F-stat p-value < 1E-8) with the
    greatest magnitude of effect from the hyper- and hypo-methylated sides, and "any", which
    selects the 100 probes (with F-stat p-value < 1E-8) with the greatest magnitude of difference
    regardless of direction of effect. Default input "auto" will use "any" for cord blood and
    "both" otherwise, in line with previous versions of this function and/or our recommendations.
    Please see the references for more details.}
  \item{cellTypes}{Which cell types, from the reference object, should be
    we use for the deconvolution? See details. }
  \item{referencePlatform}{The platform for the reference dataset; if
    the input \code{rgSet} belongs to another platform, it will be
    converted using \code{\link{convertArray}}.}
  \item{returnAll}{Should the composition table and the normalized user
    supplied data be return?}
  \item{verbose}{Should the function be verbose?}
  \item{meanPlot}{
    Whether to plots the average DNA methylation across the cell-type discrimating probes within
    the mixed and sorted samples.
}
\item{...}{Other arguments, i.e arguments passed to plots}
}
\details{
See \code{\link[minfi]{estimateCellCounts}} for more information regarding the exact details. estimateCellCounts.gds differs slightly, as it will impose the quantiles of type I and II probes onto the reference Dataset rather than normalising the two together. This is 1) More memory efficient and 2) Faster - due to not having to normalise out a very small effect the other 60 samples from the reference set will have on the remaining quantiles.

Optionally, a proportion of samples can be used to derive quantiles when there are more than 1000 samples in a dataset, this will further increase performance of the code at a cost of precision.
}
