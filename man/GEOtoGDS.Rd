\name{GEOtoGDS}
\alias{geotogds}
\alias{getpheno}
\title{
Download GEO data and pipe into a gdsfmt object
}
\description{
Uses GEOquery to download a GSE Accession into the current working directory,
only working for GSE's with raw idat files.
}
\usage{
geotogds(geo, gds, method = "wget", keepidat = F, keeptar = F, ...)
}

\arguments{
\item{geo}{
A GEO Accession number. "GSE########" or a tarball 'GSE######.tar.gz'
}
\item{gds}{
A character string specifying the name of the .gds file to write to.
}

\item{method}{
String to indicate what method to download data from GEO off. Default is 'wget'
}
\item{keepidat}{Logical, indicates whether or not to keep the raw idat files in the working directory after downloading - if FALSE, removed}
\item{keeptar}{Logical, indicates whether or not to keep the downloaded tarball in the working directory - if FALSE, removed.}
\item{...}{Additional Arguments to pass to other functions (if any)}

}

\value{
A gds.class object, which points to the newly created .gds file.
}

\author{
Tyler Gorrie-Stone
Who to contact: <tgorri@essex.ac.uk>
}

%\seealso{}
\examples{
#load example dataset
# gfile <- geotogds("GSE*******", "Nameoffile.gds")
# Will not work if gds has no idats submitted. May also fail if idats
# are not deposited in a way readily readable by readEPIC().
# closefn.gds(gfile)
}
