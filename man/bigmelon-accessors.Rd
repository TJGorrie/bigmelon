\name{bigmelon-accessors}
\alias{bigmelon-accessors}
\alias{betas,gds.class-method}
\alias{methylated,gds.class-method}
\alias{unmethylated,gds.class-method}
\alias{pvals,gds.class-method}
\alias{fData,gds.class-method}
\alias{pData,gds.class-method}
\alias{QCmethylated,gds.class-method}
\alias{QCunmethylated,gds.class-method}
\alias{QCrownames,gds.class-method}
\alias{QCmethylated}
\alias{QCunmethylated}
\alias{QCrownames}
\alias{getHistory,gds.class-method}
\alias{colnames,gds.class-method}
\alias{rownames,gds.class-method}
\alias{colnames,gdsn.class-method}
\alias{rownames,gdsn.class-method}
\alias{exprs,gds.class-method}
\alias{'[.gds.class'}
\alias{'[.gdsn.class'}
\alias{fot,gds.class-method}

\title{
Bigmelon accessors
}

\description{
Functions to access data nodes in gds.class objects.
}

\usage{

\S4method{betas}{gds.class}(object)
\S4method{methylated}{gds.class}(object)
\S4method{unmethylated}{gds.class}(object)
\S4method{pvals}{gds.class}(object)
\S4method{fData}{gds.class}(object)
\S4method{pData}{gds.class}(object)
\S4method{QCmethylated}{gds.class}(object)
\S4method{QCunmethylated}{gds.class}(object)
\S4method{QCrownames}{gds.class}(object)
\S4method{getHistory}{gds.class}(object)
\S4method{colnames}{gds.class}(x, do.NULL=TRUE, prefix=NULL)
\S4method{rownames}{gds.class}(x, do.NULL=TRUE, prefix=NULL)
\S4method{exprs}{gds.class}(object)
\S4method{fot}{gds.class}(x)
}

\arguments{

\item{object}{
A gds.class object.
}


for colnames and rownames:
\item{x}{
A gds.class object.
}

\item{do.NULL}{
logical. If 'FALSE' and names are 'NULL', names are created.
}

\item{prefix}{
prefix: for created names.
}

}

\details{
Each function returns the data stored in the corresponding node as either a
gdsn.class object or a matrix or data.frame. These are names following the
conventions of the methylumi package and perform similar functions.

Each function which returns a gdsn.class object can be indexed using
matrix-like '[' operations. With an optional name argument which optionally
allows for row and col names to be automatically appended to returned matrix.

The QC functions (returns QCdata split into separate matrices for methylated
values, unmethylated values, and probe names)

exprs returns a data.frame of beta values for all probes across all samples.

fot returns a vector corresponding to probe design (either 'I' or 'II') which is used for normalising differences between probe designs.

}

\value{
Returns specified node representing the called accessor
}

\author{
Leonard C Schalkwyk, Ayden Saffari, Tyler Gorrie-Stone
Who to contact: <tgorri@essex.ac.uk>
}

\seealso{
\code{\link{bigmelon}}, \code{\link[methylumi]{methylumi}}
}

\examples{
data(melon)
e <- es2gds(melon,'wat_melon.gds')
betas(e)
betas(e)[,]
methylated(e)[1:5, ]
unmethylated(e)[ ,1:5]
pvals(e)[1:5, 1:5]
fData(e)
pData(e)
colnames(e)
rownames(e)
exprs(e)

closefn.gds(e)
unlink("wat_melon.gds")
}
