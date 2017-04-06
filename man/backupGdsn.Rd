\name{backup.gdsn}
\alias{backup.gdsn}

\title{
Copy gds node to a backup folder within gds object
}

\description{
Quick function that will copy designated gdsn.class node within a gds object
to a 'backup'folder. If 'backup' folder does not exist, this is created. This is a wrapper to \code{\link[gdsfmt]{copyto.gdsn}} which should be used if one wishes to copy a gds node to a seperate gds file.
}

\usage{
backup.gdsn(gds = NULL, node)
}

\arguments{
\item{gds}{
If NULL, function will call \code{\link[gdsfmt]{getfolder.gdsn}} to find the
root node. Otherwise, user can specify a separate gds.class object to copy
the specified node to.
}

\item{node}{
gdsn.class object (a gds node)
}
}

\value{
gds object is modified to have a new folder 'backup' with supplied node
copied inside
}

\author{
Tyler Gorrie-Stone <tgorri@essex.ac.uk>
}

\seealso{\code{\link[gdsfmt]{copyto.gdsn}}}

\examples{
data(melon)
e <- es2gds(melon, "melon.gds")
nod <- index.gdsn(e, "betas")
backup.gdsn(gds = NULL, node = nod)

closefn.gds(e)
unlink("melon.gds")
}

\keyword{backup}
