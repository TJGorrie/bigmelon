\name{pfilter.gds}
\alias{pfilter.gds}
\alias{pfilter,gds.class-method}
\alias{pfilter,gds.class}
\title{Basic data filtering for Illumina methylation data in gds objects}
\description{
The pfilter function filters data sets based on bead count and
detection p-values.  The user can set their own thresholds or use
the default pfilter settings. This specific function will take a Genomic Data
Structure (GDS) file as input and perform pfilter similar to how
\code{\link[wateRmelon]{pfilter}} in wateRmelon is performed.
}
\usage{
\S4method{pfilter}{gds.class}(mn, perCount = NULL, pnthresh = NULL,
perc = NULL, pthresh = NULL)
}
\arguments{
\item{mn}{
a gds object OR node corresponding to methylated intensities
}
\item{perCount}{
remove sites having this percentage of samples with a beadcount <3, default = 5
}
\item{pnthresh}{
cut off for detection p-value, default= 0.05
}
\item{perc}{
remove sample having this percentage of sites with a detection p-value
greater than pnthresh, default = 1
}
\item{pthresh}{
remove sites having this percentage of sample with a detection p-value greater
than pnthresh, default = 1
}
}
%\details{}
\value{ See \code{\link[wateRmelon]{pfilter}}. If using pfilter.gds, function
will return a list of two elements. Both logical vectors with length(nrow) and length(ncol), which can be used for subsetting. Otherwise will be used to subset data directly.
}
\author{Tyler Gorrie-Stone, Original (wateRmelon) Function by Ruth Pidsley}
\seealso{\code{\link[wateRmelon]{pfilter}}}
\examples{
data(melon)
e <- es2gds(melon, "melon.gds")
pfilter(e)
closefn.gds(e)
unlink("melon.gds")
}
