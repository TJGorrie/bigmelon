\name{bumphunterEngine.gdsn}
\alias{bumphunterEngine.gdsn}
\alias{.getEstimate2}
\title{
Bumphunter using bigmelon
}
\description{
Estimate regions for which a genomic profile deviates from its baseline value. Originally implemented to detect differentially methylated genomic regions between two populations. Functions identically to \code{\link[bumphunter]{bumphunter}}. 
}
\usage{
bumphunterEngine.gdsn(mat, design, chr = NULL, pos, cluster = NULL, coef = 2, cutoff = NULL, pickCutoff = FALSE, pickCutoffQ = 0.99, maxGap = 500, nullMethod=c("permutation","bootstrap"), smooth = FALSE, smoothFunction = locfitByCluster, useWeights = FALSE, B=ncol(permutations), permutations=NULL, verbose = TRUE, ...)
}
\arguments{
  \item{mat}{A gdsn.class object (e.g betas(gfile)}
  \item{design}{Design matrix with rows representing samples and columns
    representing covariates. Regression is applied to each row of mat.}
  \item{chr}{A character vector with the chromosomes of each location.}
  \item{pos}{A numeric vector representing the chromosomal position.}
  \item{cluster}{The clusters of locations that are to be analyzed
    together. In the case of microarrays, the clusters are many times
    supplied by the manufacturer. If not available the function
    \code{\link[bumphunter]{clusterMaker}} can be used to cluster nearby locations.}
  \item{coef}{An integer denoting the column of the design matrix
    containing the covariate of interest. The hunt for bumps will be
    only be done for the estimate of this coefficient.}
  \item{cutoff}{A numeric value. Values of the estimate of the genomic
    profile above the cutoff or below the negative of the cutoff will be
    used as candidate regions. It is possible to give two separate
    values (upper and lower bounds). If one value is given, the lower
    bound is minus the value.}
  \item{pickCutoff}{Should bumphunter attempt to pick a cutoff using the
    permutation distribution?}
  \item{pickCutoffQ}{The quantile used for picking the cutoff using the
    permutation distribution.}
  \item{maxGap}{
    If cluster is not provided this maximum location gap will be used to define cluster
    via the \code{\link[bumphunter]{clusterMaker}} function.}
  \item{nullMethod}{
    Method used to generate null candidate regions, must be one of \sQuote{bootstrap} or
    \sQuote{permutation} (defaults to \sQuote{permutation}). However, if covariates in addition to the 
    outcome of interest are included in the design matrix
  (ncol(design)>2), the \sQuote{permutation} approach is not
  recommended. See vignette and original paper for more information.}
  \item{smooth}{
    A logical value. If TRUE the estimated profile will be smoothed with the
    smoother defined by \code{smoothFunction}}
  \item{smoothFunction}{
    A function to be used for smoothing the estimate of the genomic
    profile. Two functions are provided by the package: \code{loessByCluster}
    and \code{runmedByCluster}.}
  \item{useWeights}{
    A logical value. If \code{TRUE} then the standard errors of the
   point-wise estimates of the profile function will be used as weights
   in the loess smoother \code{loessByCluster}. If the
   \code{runmedByCluster} smoother is used this argument is ignored.}
 \item{B}{An integer denoting the number of resamples to use when computing
   null distributions. This defaults to 0. If \code{permutations} is
  supplied that defines the number of permutations/bootstraps and \code{B} is
  ignored.}
\item{permutations}{ is a matrix with columns providing indexes to be used to
  scramble the data and create a null distribution when
  \code{nullMethod} is set to  permutations. If the bootstrap approach is used this argument is
  ignored. If this matrix is not supplied and \code{B}>0 then
  these indexes are created using the function \code{sample}. }
\item{verbose}{
    logical value. If \code{TRUE}, it writes out some messages
    indicating progress. If \code{FALSE} nothing should be printed.
  }
  \item{\dots}{
    further arguments to be passed to the smoother functions.
}
}
\details{
This function is a direct replication of the \code{\link[bumphunter]{bumphunter}} function by Rafael A. Irizarry, Martin J. Aryee, Kasper D. Hansen, and Shan Andrews.
}
\author{
Original Function by Rafael A. Irizarry, Martin J. Aryee, Kasper D. Hansen, and Shan Andrews.
Bigmelon implementation by Tyler Gorrie-Stone
Who to contact if this all goes horribly wrong: <t.gorrie-stone@qmul.ac.uk>
}

