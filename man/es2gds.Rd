\name{es2gds}
\alias{es2gds}
\title{
Coersion method for MethyLumiSet, RGChannelSet or MethylSet to CoreArray 
Genomic Data Structure (GDS) data file
}
\description{
The es2gds function takes a MethyLumiSet, RGChannelSet or MethylSet data object
and converts it into a CoreArray Genomic Data Structure (GDS) data file 
(via the gdsfmt package), returning this as a gds.class object for use 
with bigmelon. 
}
\usage{
es2gds(m, file, qc = TRUE)
}

\arguments{
\item{m}{
A MethyLumiSet, RGChannelSet or MethylSet object
}
\item{file}{
A character string specifying the name of the .gds file to write to.
}

\item{qc}{
When set to true (default), data from control probes included.
}
}

\value{
A gds.class object, which points to the newly created .gds file.
}

\author{
Leonard C Schalkwyk, Ayden Saffari, Tyler Gorrie-Stone
Who to contact: <t.gorrie-stone@qmul.ac.uk>
}

\seealso{
\code{\link{app2gds}}, \code{\link{iadd}}.

}
\examples{
#load example dataset
data(melon)
#convert to gds
e <- es2gds(melon,'melon.gds')
closefn.gds(e)
unlink('melon.gds')
}

%% Add one or more standard keywords, see file 'KEYWORDS' in the
%% R documentation directory.
%%\keyword{ ~kwd1 }
%%\keyword{ ~kwd2 }
