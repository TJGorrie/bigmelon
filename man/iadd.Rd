\name{iadd}
\alias{iadd}
\alias{iadd2}

\title{
Add data from IDAT files for a single barcode to a gds file

}
\description{
Add data from IDAT files for a single barcode to a gds file.
or 
Add data from many IDAT files from a single directory to a gds file.
}
\usage{
iadd(bar, gds,...)
iadd2(path, gds, chunksize = NULL, ...)
}

\arguments{
\item{bar}{
The barcode for an IDAT file Or the file path of the file containing red or 
green channel intensities for that barcode (this will automatically 
locate and import both files regardless of which one you provide)
}
\item{path}{
The file path where (multiple) IDAT files exist. Iadd2 will process every idat 
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
\item{...}{
Additional Arguments passed to \code{\link[wateRmelon]{methylumIDATepic}}.
}
}

\value{
returns a gds.class object, which points to the appended .gds file.
}

\author{
Tyler Gorrie-Stone, Leonard C Schalkwyk, Ayden Saffari.
Who to contact: <tgorri@essex.ac.uk>
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
