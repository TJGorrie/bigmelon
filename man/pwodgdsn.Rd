\name{pwod.gdsn}
\alias{pwod.gdsn}
\alias{pwod, gds.class}
\alias{pwod, gdsn.class}

\title{
Probe-Wise Outlier Detection
}
\description{
'P'robe-'W'ise 'O'utlier 'D'etection via interquartile ranges
}
\usage{
pwod.gdsn(node, mul = 4)
}
\arguments{
\item{node}{
gdsn.class object containing array to be filtered
}
\item{mul}{
Number of interquartile ranges used to determine outlying
probes. Default is 4 to ensure only very obvious outliers are
removed.
}
}
\details{
Detects outlying probes across arrays in methylumi and minfi
objects. Outliers are probable low MAF/SNP heterozygotes.
}
\value{
Nothing returned. Supplied gds object will have new node with outlier
probes coerced to NAs.
}

\author{
Tyler Gorrie-Stone
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
