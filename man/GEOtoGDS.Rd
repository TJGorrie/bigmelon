\name{GEOtoGDS}
\alias{geotogds}
\alias{getpheno}
\title{
Download data from GEO and convert it into a gdsfmt object
}
\description{
Uses the GEOquery R package to download a GSE Accession into the current working directory. This will only work for GSE's that have raw idat files associated with them.
}
\usage{
geotogds(geo, gds, method = "wget", keepidat = F, keeptar = F, ...)
}

\arguments{
\item{geo}{
Either a GEO accession number ('GSE########') or a previously downloaded tarball 'GSE######.tar.gz'
}
\item{gds}{
A character string that specifies the path and name of the .gds file you want to write to.
}
\item{method}{
Character value to indicate which method should be used to download data from GEO. Default is 'wget'
}
\item{keepidat}{Logical, indicate whether or not raw idat files in the working directon should be removed after parsing, if FALSE: idat files will be removed.}
\item{keeptar}{Logical, indicate whether or not the downloaded tarball should be removed after parsing, if FALSE: the tarball will be removed.}
\item{...}{
Additional Arguments to pass to other functions (if any)}
}
\value{
geotogds will return a gds.class object that will point towards a the newly created .gds file with majority of downloaded contents inside.
}

\author{
Tyler Gorrie-Stone
Who to contact: <t.gorrie-stone@qmul.ac.uk>
}

%\seealso{}
\examples{
#load example dataset
# gfile <- geotogds("GSE*******", "Nameoffile.gds")
# Will not work if gds has no idats submitted. May also fail if idats
# are not deposited in a way readily readable by readEPIC().
# closefn.gds(gfile)
}
