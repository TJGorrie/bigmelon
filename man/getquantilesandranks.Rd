\name{getquantilesandranks}
\alias{getquantilesandranks}
\alias{.sortvector}
\title{
Compute the quantiles and ranks for a given gdsn.node
}
\description{
Used inside \code{\link{dasenrank}} to generate the quantiles for both type 'I' and type 'II' probes to normalise DNA methylation data using bigmelon.
}
\usage{
getquantilesandranks(gds, node, onetwo, rank.node = NULL, perc = 1)
}
\arguments{
\item{gds}{
A gds.class object
}
\item{node}{
A gdsn.class object, or a character string that refers to a node within supplied gds.
}
\item{onetwo}{
gdsn.class object OR character string that refers to location in gds that contains information relating to probe design OR vector of length equal to the number of rows in the array that contains 'I' and 'II' in accordance to Illumina HumanMethylation micro-array design. This can be obtained with fot(gds)
}
\item{rank.node}{Default = NULL. If supplied with character string, function will calculate the ranks of given node and store them in gds. Additionally, the computed quantiles will now instead be attributed to rank.node which can be accessed with \code{\link{get.attr.gdsn}}
}
\item{perc}{
A number between 0 and 1 that relates to the given proportion of columns that are used to normalise the data. Default is set to 1, but in cases where there many of samples to normalise this number can be reduced to increase speed of code.
}
}
\details{
Used in \code{\link{dasenrank}}, can be used externally for testing purposes.
}
\value{
If rank.node is NULL. A list containing quantiles, intervals and supplied probe design will be returned.
If rank.node was supplied, nothing will be returned. Instead a new node will be created in given gds that has the otherwise returned list attached as an attribute. Which can be accessed with \code{\link{get.attr.gdsn}}
}

\author{
Tyler Gorrie-Stone
Who to contact: <t.gorrie-stone@qmul.ac.uk>
}
%\seealso{
%}
\examples{
data(melon)
e <- es2gds(melon, "melon.gds")
output <- getquantilesandranks(gds = e, 'methylated', onetwo = fot(e), perc = 1, rank.node = NULL)
# with-out put.
#getquantilesandranks(gds = e, 'methylated', onetwo = fot(e), perc = 1, rank.node = 'mnsrank')
closefn.gds(e)
unlink("melon.gds")
}
