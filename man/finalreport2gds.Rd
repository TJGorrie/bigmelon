\name{finalreport2gds}
\alias{finalreport2gds}
\title{
Read finalreport files and convert to genomic data structure files
}
\description{
Function to easily load Illumina methylation data into a genomic data 
structure (GDS) file.
}
\usage{
finalreport2gds(finalreport, gds, ...)
}
\arguments{
\item{finalreport}{
A filename of the text file exported from GenomeStudio
}
\item{gds}{
The filename for the gds file to be created
}
\item{\dots}{
Additional arguments passed to \code{\link[methylumi]{methylumiR}}
}
}
\details{
Creates a .gds file.
}
\value{
A gds.class object
}

\author{
Tyler Gorrie-Stone
}

\examples{
finalreport <- "finalreport.txt"
\dontrun{finalreport2gds(finalreport, gds="finalreport.gds")}
}

