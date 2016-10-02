\name{app2gds}
\alias{app2gds}
\title{
Append a MethyLumiSet object to a CoreArray Genomic Data Structure (GDS) data 
file
}
\description{
This function will append a MethyLumiSet data object to a CoreArray Genomic 
Data Structure(GDS) data file, and return as a gds.class object.
}
\usage{
app2gds(m, bmln)
}

\arguments{
\item{m}{
The MethyLumiSet object to be appended to the CoreArray Genomic Data Structure 
(GDS) data file
}
\item{bmln}{

Either: 
A gds.class object

Or:
A character string specifying the name of an existing .gds file to write to.

Or:
A character string specifying the name of a new .gds file to write to
}
}

\value{
A gds.class object, which points to the appended .gds file.
}

\author{
Leonard C Schalkwyk, Ayden Saffari, Tyler Gorrie-Stone
Who to contact: <tgorri@essex.ac.uk>
}

\seealso{
\code{\link{es2gds}}, \code{\link{iadd}}.
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
