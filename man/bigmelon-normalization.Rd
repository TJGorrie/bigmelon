\name{bigmelon-normalization}
\alias{dasen.gds}
\alias{dasen,gds.class-method}
\alias{danen.gds}
\alias{danen,gds.class-method}
\alias{naten.gds}
\alias{naten,gds.class-method}
\alias{daten1.gds}
\alias{daten1,gds.class-method}
\alias{daten2.gds}
\alias{daten2,gds.class-method}
\alias{nasen.gds}
\alias{nasen,gds.class-method}
\alias{nanet.gds}
\alias{nanet,gds.class-method}
\alias{danet.gds}
\alias{danet,gds.class-method}
\alias{nanes.gds}
\alias{nanes,gds.class-method}
\alias{danes.gds}
\alias{danes,gds.class-method}
\alias{qn.gdsn}
\alias{db.gdsn}
\alias{dfsfit.gdsn}
\alias{betaqn,gds.class-method}

\title{
Bigmelon Quantile Normalization methods.
}

\description{
Set of functions that are used to perform quantile normalization methods on gds.class objects.
}

\usage{
\S4method{dasen}{gds.class}(mns, fudge = 100, ret2 = FALSE, node="betas",...)
dasen.gds(gds, node, mns, uns, onetwo, roco, fudge, ret2)
qn.gdsn(gds, target, newnode)
db.gdsn(gds, mns, uns)
dfsfit.gdsn(gds, targetnode, newnode, roco, onetwo)
}

\arguments{
\item{gds}{ A gds.class object
}
\item{node}{ The "name" of desired output \code{\link[gdsfmt]{gdsn.class}} node
}
\item{mns}{ The \code{gdsn.class} node that corresponds to "methylated"
intensities.
}
\item{uns}{ The \code{gdsn.class} node that corresponds to "unmethylated"
intensities.
}
\item{onetwo}{ The \code{gdsn.class} node that corresponds to probe designs
(in reference to 450k and EPIC arrays)
OR
character string pointing to location of \code{gdsn.class} node. e.g.
"fData/DESIGN"
OR
vector containing probe design types of length > 1.
}
\item{roco}{ This allows a background gradient model to be fit.
This is split from data column names by default.  roco=NULL
disables model fitting (and speeds up processing), otherwise
roco can be supplied as a character vector of strings like
'R01C01' (only 3rd and 6th characters used).
}
\item{fudge}{ The value added to total intensity to prevent denominactors close
to zero when calculation betas. default = 100
}
\item{ret2}{ if TRUE, appends the newly calculated methylated and
unmethylated intensities to original gds (as specified in gds arguement).
Will overwrite the raw intensities.
}
\item{target}{ Target \code{gdsn.class} node to perform normalization on.
If using "*****.gds" method you do not need to specify this.
}
\item{targetnode}{ Target \code{gdsn.class} node to perform normalization
on. If using "*****.gds" method you do not need to specify this.
}
\item{newnode}{ "name" of desired output \code{gdsn.class} node.
If using "*****.gds" method you do not need to specify this.
}
\item{...}{Additional args such as roco or onetwo.}
}

\details{
Each function performs a normalization method described within the wateRmelon
package. Functions: \code{qn.gdsn}, \code{design.qn.gdsn}, \code{db.gdsn} and
\code{dfsfit.gdsn} are described to allow users to create their own custom
normalization methods. Otherwise calling \code{dasen} or \code{dasen.gds}
e.t.c will perform the necessary operations for quantile normalization.

Each 'named' normalization method will write a temporary gds object called
"temp.gds" into the current working directory and it is removed
when normalization is complete. Current methods supplied by default arguments
will replace the raw betas with normalized betas, but leave the methylated and unmethylated intensities unprocessed.
}

\value{
Normalization methods return nothing but will affect the gds file and
replace/add nodes given to the function.
}

\author{
Tyler J Gorrie-Stone <t.gorrie-stone@qmul.ac.uk>
}
\seealso{
\code{\link{wateRmelon}}, \code{\link[wateRmelon]{dasen}}
}
\examples{
data(melon)
e <- es2gds(melon,'wat_melon.gds')
dasen(e)
closefn.gds(e) # Close gds object
unlink('wat_melon.gds') # Delete Temp file
}
