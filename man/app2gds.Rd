\name{app2gds}
\alias{app2gds}
\title{
Append a MethyLumiSet object to a preexisting gds file
}
\description{
This function will append a MethyLumiSet object to a gds file and return a gds.class object.
}
\usage{
app2gds(m, bmln)
}

\arguments{
\item{m}{
The MethyLumiSet object to be appended to the gds file, with the same number of rows as the gds file.
}
\item{bmln}{

Either:
A gds.class object

Or:
A character string specifying the filepath of an existing .gds file to write to.

Or:
A character string specifying the file path of a new .gds file to write to
}
}

\details{
Currently this function only takes a MethyLumiSet object as the only type of eligible input. This function will also produce unexpected results if the number of rows of the new object does not match the existing .gds file. Hopefully the function will noisily fail if this is the case however to prevent any errors from occuring it is recommended that raw idat files are read in using \code{\link[wateRmelon]{readEPIC}} or appended with \code{\link{iadd}} or \code{\link{iadd2}} to ensure that all rows are of the same length and have the same annotation.
}

\value{
A gds.class object pointed towards the newly appended .gds file.
}

\author{
Leonard C Schalkwyk, Ayden Saffari, Tyler Gorrie-Stone
Who to contact: <t.gorrie-stone@qmul.ac.uk>
}

\seealso{
\code{\link{es2gds}}, \code{\link{iadd}, \code{\link{iadd2}}.
}

\examples{
#load example dataset
data(melon)

#split data into halves
melon_1 <- melon[,1:6]
melon_2 <- melon[,7:12]

#convert first half to gds
e <- es2gds(melon_1,'1_half_melon.gds')

#append second half to existing gds file
f <- app2gds(melon_2,e)

closefn.gds(e)
unlink("1_half_melon.gds")
}

%% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }
