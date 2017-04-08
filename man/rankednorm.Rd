\name{dasenrank}
\alias{dasenrank}
\alias{computebeta.gds}
\title{
Dasen Quantile Normalization by storing ranks
}
\description{
Performs 'dasen' normalization for .gds format objects while storing
the ranked methylated/unmethylated intensities until they are needed down stream.

This will eliminate the secondary re-sorting required by quantile normalisation as it is perfomed down-stream - either by using \code{computebeta.gds} or manually with '['
}
\usage{
dasenrank(gds, mns, uns, onetwo, roco, calcbeta = NULL, perc = 1)
computebeta.gds(gds, new.node, mns, uns, fudge)
}
\arguments{
\item{gds}{
gds.class object which contains methylated and unmethylated intensities. The function will write two(four) nodes to this object called 'mnsrank' and 'unsrank' which contain the ranks of the given nodes.
}
\item{mns}{
gdsn.class object OR character string that refers to location in gds that relates to the (raw) methylated intensities.
}
\item{uns}{
gdsn.class object OR character string that refers to location in gds that relates to the (raw) unmethylated intensities.
}
\item{onetwo}{
gdsn.class object OR character string that refers to location in gds that contains information relating to probe design OR vector of length equal to the number of rows in the array that contains 'I' and 'II' in accordance to Illumina HumanMethylation micro-array design.
}
\item{roco}{
Sentrix (R0#C0#) position of all samples.
}
\item{calcbeta}{
Default = NULL, if supplied with a string, a new gdsn.node will be made with supplied string, which will contain the calculated betas.
}
\item{perc}{
A number between 0 and 1 that relates to the given proportion of columns that are used to normalise the data. Default is set to 1, but incase there are lots of samples to normalise this number can be reduce to increase speed of code.
}
\item{new.node}{
Character string depicting name of new betas node in given gds object.
}
\item{fudge}{ Arbitrary value to offset low intensities}
}
\details{
calcbeta is a known bottle-neck for this code! Also function is highly experimental.
}
\value{
Nothing is returned to the R environment, however the supplied gds will have 4 or 5 gdsn.nodes added. These are: 'mnsrank', 'unsrank', 'isnamnsrank' (hidden), 'isnaunsrank'(hidden) and calcbeta if supplied. 'mnsrank' and 'unsrank' have been given some attributes - which contain the calculated quantiles from \code{getquantilesandranks}.
}

\author{
Tyler J. Gorrie-Stone - tgorri at essex.ac.uk
}
%\seealso{
%}
\examples{
data(melon)
e <- es2gds(melon, "melon.gds")
#dasenrank(gds = e)
closefn.gds(e)
unlink("melon.gds")
}
