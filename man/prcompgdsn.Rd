\name{prcomp.gdsn.class}
\alias{prcomp.gdsn.class}
\alias{prcomp.gdsn}
\alias{prcomp, gdsn.class}
\alias{prcomp, gdsn.class-method}

\title{
Principal Component Analysis for high-dimensional data
}
\description{
Performs principal components analysis on the given
gds object and returns the results as an object of class "prcomp".
}
\usage{
\method{prcomp}{gdsn.class}(x, center = FALSE, scale. = FALSE, tol = NULL, 
rank. = NULL, retx = FALSE, perc = 0.01, npcs = NULL, ...)
}

\arguments{
\item{x}{
GDS Node
}
\item{center}{
Logical value indicating whether variables should be shifted to be zero 
centered.
}
\item{scale.}{
Logical value indicating whether the variables should be scaled to have unit 
variance
}
\item{tol}{
a value indicating the magnitude below which components
should be omitted.
}
\item{rank.}{
Old. (Still functional) Number of principal components to be returned
}
\item{retx}{
a logical value indicating whether the rotated variables
should be returned.
}
\item{perc}{
Percentage of data to be used. 
}
\item{npcs}{
Number of principal components to be returned
}
\item{\dots}{
arguments passed to or from other methods. If "x" is a
formula one might specify "scale." or "tol".
}
}
\details{
The calculation is done by a singular value decomposition of the
(centered and possibly scaled) data matrix, not by using "eigen"
on the covariance matrix. This is generally the preferred method
for numerical accuracy. The "print" method for these objects
prints the results in a nice format and the "plot" method produces
a scree plot.
}
\value{
An object of prcomp class
}

\examples{
data(melon)
e <- es2gds(melon, "melon.gds")
prcomp(betas(e))
closefn.gds(e)
unlink("melon.gds")
}

