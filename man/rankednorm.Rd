\name{dasenrank}
\alias{dasenrank}
\alias{quickquan}
\alias{computebetas}
\title{
Quantile Normalization storing ranks
}
\description{
Performs 'dasen' normalization for .gds format objects by storing the ranked 
values. This elimates the secondary reading of data until data is needed to be
read down-stream.
}
\usage{
dasenrank(gds, mns, uns, onetwo, roco, calcbeta = NULL, ...)
}
\arguments{
\item{gds}{
gds.class object containing mns and uns.
}
\item{mns}{
gdsn.class node that contains information relating to raw methylated intensity
}
\item{uns}{
gdsn.class node that contains information relating to raw unmethylated intensity
}
\item{onetwo}{
gdsn.class node that contains information relating to probe design
}
\item{roco}{
Sentrix (R0#C0#) position of all samples
}
\item{calcbeta}{
Default = NULL, if supplied with a string, a new gdsn.node will be made with supplied string.
}
\item{...}{
Further Args
}
}
\details{
May be buggy
}
\value{
Nothing returned. Supplied gds object will have new node with outlier
probes coerced to NAs
}

\author{
Tyler J. Gorrie-Stone - tgorri at essex.ac.uk
}
\seealso{
\code{\link[wateRmelon]{pwod}}
}
\examples{
data(melon)
e <- es2gds(melon, "melon.gds")
#dasenrank(gds = e)
closefn.gds(e)
unlink("melon.gds")
}

