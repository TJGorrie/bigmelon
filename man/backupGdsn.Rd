\name{backup.gdsn}
\alias{backup.gdsn}

\title{
Copy gds node to a backup folder within gds object
}

\description{
This function will copy a designated gdsn.class object stored inside a gds object to a backup folder (aptly named backup). If the backup folder does not exist, this will be created. This is a wrapper to \code{\link[gdsfmt]{copyto.gdsn}} which should be used if one wishes to copy a gds node to a seperate gds file.
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
gdsn.class object (a gds node) which can be specified using \code{\link[gdsfmt]{index.gdsn}}
}
}

\value{
A gdsobject with an additional folder called backup with supplied node within.
}

\author{
Tyler Gorrie-Stone <t.gorrie-stone@qmul.ac.uk>
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
