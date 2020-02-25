\name{gds2mlumi}
\alias{gds2mlumi}
\alias{gds2mset}

\title{
Convert Genomic Data Structure file to Methylumiset or Methylset object.
}
\description{
Convert a Genomic Data Structure object back into a methylumi object,
with subsetting features.
}
\usage{
gds2mlumi(gds, i, j)
gds2mset(gds, i, j, anno)
}

\arguments{
\item{gds}{
a gds object
}
\item{i}{
Index of rows
}
\item{j}{
Index of Columns
}
\item{anno}{
If NULL, function will attempt to guess the annotation to be used. 
Otherwise can be specified with either "27k", "450k", "epic" or "unknown".
}
}
\value{
A methylumi object
}

\author{
Tyler Gorrie-Stone
Who to contact: <t.gorrie-stone@qmul.ac.uk>
}
\examples{
data(melon)
e <- es2gds(melon, "melon.gds")
gds2mlumi(e)
closefn.gds(e)
unlink("melon.gds")
}

