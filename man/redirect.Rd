\name{redirect.gds}
\alias{redirect.gds}

\title{
Change the location of file 'paths' for row and column names.
}

\description{
Quickly change contents of gdsn.class node "paths".
}

\usage{
redirect.gds(gds, rownames, colnames)
}

\arguments{
\item{gds}{ \code{gds.class} object containing node named "paths".
}
\item{rownames}{ Character string that points to named part of supplied gds
that corresponds to rownames.
e.g. "fData/Target_ID". Default = "fData/Probe_ID"
}
\item{colnames}{ Character string that points to names part of supplied gds 
that corresponds to colnames.
e.g. "pData/Sample_ID". Default = "pData/barcode"
}
}

\value{
Changes the gdsn.class object named "paths" to supplied rownames and colnames 
within supplied gds.class object.
}

\details{
Function is important within many functions and can lead to errors if row and 
column names are not correctly specified. By default, es2gds can recognize 
whether a methylumiset object was read in through readEPIC or methylumiR 
and will set the row and col names paths correctly.
Will fail noisily if given a pathway that does not exist.
}

\author{
Tyler J. Gorrie-Stone
Who to contact: <tgorri@essex.ac.uk>
}

\seealso{
\code{\link[gdsfmt]{add.gdsn}}, \code{\link{app2gds}}, \code{\link{es2gds}}
}

\examples{
data(melon)
e <- es2gds(melon, "melon.gds") # Create gds object
redirect.gds(e, rownames = "fData/TargetID", colnames = "pData/sampleID")
# Deleting Temp files
closefn.gds(e)
unlink("melon.gds")

}
