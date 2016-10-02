\name{combo.gds}
\alias{combo.gds}

\title{
Combine two different gds objects.
}

\description{
Combines the shared gdsn.class nodes between two gds objects depending on 
primary gds.object dimensions.
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
Will crudely combine shared nodes between primary and secondary based on the 
dimensions / rownames of the primary node. NAs will be coerced where probes 
are missing from secondary gds.

Will only look for nodes with the names "betas", "methylated", "unmethylated", 
"pvals" and "NBeads".
}

\note{
Will lose information relating to "pData". Therefore we recommend compiling 
separate pData object manually.
}

\author{
Tyler Gorrie-Stone <tgorri@essex.ac.uk>
}

\value{
a new gds object that has both files within it
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
