\name{lifetab2}
\alias{lifetab2}
\alias{plot.lifetab2}
\alias{lines.lifetab2}
\alias{points.lifetab2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Formula wrapper for \code{\link{lifetab}} from the
\code{KMsurv} package.
}
\description{
Calculate a life table using the actuarial method using a formula and a
data-frame with optional breaks. 
}
\usage{
lifetab2(formula, data, subset, breaks = NULL)
\method{plot}{lifetab2}(x, y=NULL, type="s", ...)
\method{lines}{lifetab2}(x, y=NULL, type="s", ...)
\method{points}{lifetab2}(x, y=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
    formula with the left-hand side being a \code{Surv} object, including a time and event indicator, and the right-hand side indicated stratification.
}
  \item{data}{
    optional \code{data.frame} for the \code{Surv} object. If this is not provided, then the parent frame
    is used for the \code{Surv} object.
}
  \item{subset}{
    optional \code{subset} statement
}
  \item{breaks}{
    optional numeric vector of breaks. If this is not provided, then the
    unique time values from the \code{Surv} object are used together
    with \code{Inf}.
  }
  \item{x}{\code{lifetab2} object}
  \item{y}{unused argument (part of the generic function)}
  \item{type}{line type}
  \item{\dots}{other arguments}
}
\details{
See \code{lifetab} for details. This wrapper is meant to make
life easier.

A copy of the \code{lifetab} function has been included in
the \pkg{biostat3} package to reduce dependencies.
}
\value{
  A \code{data.frame} as per \code{\link{lifetab}}.
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{
Mark Clements for the wrapper.
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
## we can use unique transformed times (colon_sample)
lifetab2(Surv(floor(surv_yy),status=="Dead: cancer")~1, colon_sample)

## we can also use the breaks argument (colon)
lifetab2(Surv(surv_yy,status=="Dead: cancer")~1, colon, breaks=0:10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }% use one of  RShowDoc("KEYWORDS")

