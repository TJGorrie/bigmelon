\name{redirect.gds}
\alias{redirect.gds}

\title{
Change the location of the paths for row and column names in a gds file.
}

\description{
Change how a gds object ascribes row and column names by changing the paths node of a gds file.
Although a bit laborious, the row and column information is not preserved inside a node which contains data. Thusly
during accession the row and column names are attributed after data has been accessed.
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
Changes the gdsn.class node named "paths" to supplied rownames and colnames 
within supplied gds.class object.
}

\details{
This function is particularly important within many functions inside bigmelon and
may lead to downstream errors if the row and column names are not correctly specified. 
If data is read in through es2gds the path nodex should be correctly set up and all downstream
analysis will carry out as normal.
Will fail noisily if given a pathway that does not exist.
}

\author{
Tyler J. Gorrie-Stone
Who to contact: <t.gorrie-stone@qmul.ac.uk>
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
