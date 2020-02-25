\name{combo.gds}
\alias{combo.gds}

\title{
Combine two different gds objects together
}

\description{
Function will attempt to combine together the shared gdsn.class nodes between two gds object 
depending on the dimensions of the primary gds.class object.
}

\usage{
combo.gds(file, primary, secondary)
}

\arguments{
\item{file}{
Name of the new gds file to be created.
}
\item{primary}{
A gds.class object.
}
\item{secondary}{
A gds.class object.
}
}

\details{
--EXPERIMENTAL--
Will crudely combine shared nodes between primary and secondary based on the
dimensions / rownames of the primary node. NAs will be coerced where probes
are missing from secondary gds.

Currently will only look for nodes with the names "betas", "methylated", "unmethylated",
"pvals" and "NBeads".
}

\note{
Will lose information relating to "pData". Therefore we recommend compiling
separate pData object manually and adding combined pData post-function
}

\author{
Tyler Gorrie-Stone <t.gorrie-stone@qmul.ac.uk>
}

\value{
Returns (and creates) as new gds file in the specified location with the combination
of two gds objects together.
}

\examples{
data(melon)
a <- es2gds(melon[,1:6], "primary.gds")
b <- es2gds(melon[,7:12], "secondary.gds")

ab <- combo.gds("combo.gds", primary = a, secondary = b)

closefn.gds(a)
unlink("primary.gds")
closefn.gds(b)
unlink("secondary.gds")
closefn.gds(ab)
unlink("combo.gds")
}

\keyword{combo}
