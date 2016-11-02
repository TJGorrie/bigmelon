\name{prcomp.gds.class}
\alias{prcomp.gds.class}
\alias{prcomp.gds}
\alias{prcomp, gds.class}
\alias{prcomp, gds.class-method}

\title{
Principal Component Analysis for high-dimensional data
}
\description{
Performs principal components analysis on the given
gds object and returns the results as an object of class "prcomp".
}
\usage{
\method{prcomp}{gds.class}(x, node.name, center = FALSE, scale. = FALSE, 
rank. = NULL, retx = FALSE, tol = NULL, perc = 0.01,
npcs = NULL, parallel = NULL, method = c('quick', 'sorted'), verbose = FALSE, ...)
}

\arguments{
\item{x}{
A gds.class object.
}
\item{node.name}{
Name of the gdsn.class node to learn principal components from
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
\item{parallel}{
Can supply either a cluster object (made from makeCluster) or a integer depicting number of cores to be used.
This is only used if method="sorted".
}
\item{method}{
Indicates what method to use out of "quick" and "sorted". "quick" stochastically selects number of rows according to perc. And the supplies them to svd.
"sorted" determines the interquartile range for each row then selects the top percentage (according to perc) and supplies selected rows to svd.
}
\item{verbose}{
A logical value indicating whether message outputs are displayed.
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
prcomp(e, node.name="betas", perc=0.01, method='quick')
closefn.gds(e)
unlink("melon.gds")
}

