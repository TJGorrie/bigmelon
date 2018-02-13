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
\alias{agep, gds.class-method}
\alias{agep, gdsn.class-method}
\alias{outlyx, gds.class-method}
\alias{outlyx, gdsn.class-method}
\alias{bscon, gds.class-method}
\alias{bscon, gdsn.class-method}
\alias{qual,gdsn.class,gdsn.class-method}

\title{
Functions ported from wateRmelon
}
\description{
methods for extraenous functions from wateRmelon see respective manual pages.
}
\usage{
\S4method{qual}{gdsn.class,gdsn.class}(norm, raw)
\S4method{agep}{gds.class}(betas, coeff)
\S4method{agep}{gdsn.class}(betas, coeff)
}
\arguments{
\item{betas}{A gds.class or gdsn.class object. Ideally containing normalised betas values
}
\item{coeff}{ In wateRmelon, can be set to NULL. In bigmelon must be defined. Run data(coef) before using agep.
}

\item{norm}{
normalized node (gdsn.class)
}
\item{raw}{
raw node (gdsn.class)
}
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
