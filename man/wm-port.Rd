\name{wm-port}
\alias{wm-port}
\alias{dmrse,gds.class-method}
\alias{dmrse,gdsn.class-method}
\alias{dmrse_col,gds.class-method}
\alias{dmrse_col,gdsn.class-method}
\alias{dmrse_row,gds.class-method}
\alias{dmrse_row,gdsn.class-method}
\alias{genki,gds.class-method}
\alias{genki,gdsn.class-method}
\alias{seabi,gds.class-method}
\alias{qual,gdsn.class,gdsn.class-method}
\alias{predictSex,gdsn.class-method}
\alias{predictSex,gds.class-method}
\alias{predictSex}

\title{
Functions ported from wateRmelon
}
\description{
methods for extraenous functions from wateRmelon see respective manual pages.
}
\usage{
\S4method{qual}{gdsn.class,gdsn.class}(norm, raw)
\S4method{predictSex}{gds.class}(x, x.probes, pc=2, plot=FALSE, irlba=TRUE, center=FALSE, scale.=FALSE)
\S4method{predictSex}{gdsn.class}(x, x.probes, pc=2, plot=FALSE, irlba=TRUE, center=FALSE, scale.=FALSE)
}
\arguments{
\item{norm}{
normalized node (gdsn.class)
}
\item{raw}{
raw node (gdsn.class)
}
\item{x}{gdsclass object or node corresponding to betas}
\item{x.probes}{logical or numeric vector containing indicies of X chromosome probes}
\item{pc}{The principal component to guess sex across (usually the 2nd one in most cases)}
\item{plot}{Logical, indicated whether or not to plot the prediction}
\item{irlba}{Logical, indicate whether or not to use the faster method to generate principal components}
\item{center}{Logical, indicate whether or not to center data around 0}
\item{scale.}{Logical, indicate whether or not to scale data prior to prcomp}
}

\details{
Filler Text
}
\value{
Returns expected output of functions from wateRmelon
}

\seealso{
\code{\link{wateRmelon}}
}

%%\examples{
%%}

\keyword{wm-port}
