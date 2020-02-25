\name{pwod.gdsn}
\alias{pwod.gdsn}
\alias{pwod, gds.class}
\alias{pwod, gdsn.class}

\title{
Probe-Wise Outlier Detection for DNA methylation data.
}
\description{
Function performed outlier detection for each probe (row) using Tukey's Interquartile Range method.
}
\usage{
pwod.gdsn(node, mul = 4)
}
\arguments{
\item{node}{
gdsn.class node that contains the data matrix to be filtered
}
\item{mul}{
The number of interquartile ranges used to determine outlying
probes. Default is 4 to ensure only very obvious outliers are
removed.
}
}
\details{
Detects outlying probes across arrays in methylumi and minfi
objects.
}
\value{
Nothing is returned. However the supplied gds object (of-which the node is a child of) will have a new node
with NAs interdispersed where outliers are found.
}

\author{
Tyler Gorrie-Stone
Who to contact: <t.gorrie-stone@qmul.ac.uk>
}

\seealso{
\code{\link[wateRmelon]{pwod}}
}
\examples{
data(melon)
e <- es2gds(melon, "melon.gds")
pwod(e)
closefn.gds(e)
unlink("melon.gds")
}
