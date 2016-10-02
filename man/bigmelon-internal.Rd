\name{bigmelon-internal}
\alias{bigmelon-internal}
\alias{newgds}
\alias{handle}
\alias{findgdsobj}
\alias{bmrownames}
\alias{bmrownames,gds.class-method}
\alias{subSet}
\alias{subSet,gds.class-method}
\alias{dasen2}
\alias{dfsfit.gdsn.class.1}
\alias{normalizeQuantiles.gdsn.class}
\alias{zot}
\alias{getTrait}
\alias{[.gds.class}
\alias{[.gdsn.class}


\title{
Internal bigmelon functions

}
\description{
Internal bigmelon functions, not intended for end user.
}
\usage{

newgds(file)
handle(gds)
findgdsobj(gds)
subSet(x, i, j, ..., drop = FALSE)

}

\arguments{

\item{file}{
A character string specifying the name of the .gds file to write to.
}

\item{gds}{
Either a gds.class object, or a character string specifying an existing 
.gds file.
}

\item{x}{
A gds.class object
}

\item{i}{
rows (probes) to select for subsetting
}

\item{j}{
columns (samples) to select for subsetting
}

\item{drop}{
passed on to '[' indexing operator.
}

}

\details{
newgds creates a new .gds file stub with the given name.
handle is used by app2gds to return a file handle for the given object.
findgdsobj is used by handle to search the workspace for a gds.class object 
linked to the file name specified.
subSet is used by internal functions to select a subset of rows (probes) and 
columns (samples) from a .gds file (overwriting the existing).

}

\value{
Nothing.
}

\author{
Leonard C Schalkwyk, Ayden Saffari, Tyler Gorrie-Stone
Who to contact: <tgorri@essex.ac.uk>
}

\seealso{
\code{\link{bigmelon}}
}

\examples{
data(melon)
gfile <- es2gds(melon, "melon.gds")
betas(gfile)[1:5,1:3]
closefn.gds(gfile)
unlink("melon.gds")
}

\keyword{internal}
