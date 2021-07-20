\name{iadd}
\alias{iadd}
\alias{iadd2}
\alias{idats2gds}

\title{
Add data from Illumina IDAT files to a gds file.
}
\description{
iadd will add data from multiple, specified, idat files providing '/path/to/barcode' is valid path
to a specified gds file. Barcode here implies the first part of the idat file name i.e without '_(Red|Grn).idat'

iadd2 will add data from all idat files that are stored within a single directory to a gds file.
}
\usage{
iadd(bar, gds, n = TRUE, ...)
iadd2(path, gds, chunksize = NULL, ...)
idat2gds(barcodes, gds, n=TRUE, force = FALSE, ...)
}

\arguments{
\item{bar}{
The barcode for an IDAT file OR the file path of the file containing red or
green channel intensities for that barcode (this will automatically
locate and import both files regardless of which one you provide)
}
\item{path}{
The file path where (multiple) IDAT files exist. iadd2 will process every idat
within the specified directory.
}
\item{gds}{

Either:
A gds.class object

Or:
A character string specifying the name of an existing .gds file to write to.

Or:
A character string specifying the name of a new .gds file to write to

}
\item{chunksize}{
If NULL, iadd2 will read in all barcodes in one go. Or if supplied with
a numeric value, iadd2 will read in that number of idat files in batches
}

\item{n}{
Logical, whether or not bead-counts are extracted from idat files.
}

\item{Force}{
Logical, whether or not rownames from the first idat file are applied to all idat files. Useful in combining together idat files of differing lengths.
}

\item{...}{
Additional Arguments passed to wateRmelons methylumIDATepic.
}
}

\value{
returns a gds.class object, which points to the appended .gds file.
}

\author{
Tyler Gorrie-Stone, Leonard C Schalkwyk, Ayden Saffari.
Who to contact: <tyler.gorrie-stone@diamond.ac.uk>
}

\seealso{
\code{\link{es2gds}}, \code{\link{app2gds}}.
}

\examples{
if(require('minfiData')){
bd <- system.file('extdata', package='minfiData')
gfile <- iadd2(file.path(bd, '5723646052'), gds = 'melon.gds')
closefn.gds(gfile)
unlink('melon.gds')
}
}

%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
