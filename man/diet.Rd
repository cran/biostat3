\name{diet}
\alias{diet}
\docType{data}
\title{
Diet data set
}
\description{
Dataset on a cohort investigating the incidence of coronary heart
disease (CHD)
}
\usage{data("diet")}
\format{
  A data frame with 337 observations on the following 15 variables.
  \describe{
    \item{\code{id}}{a numeric vector individual identifier}
    \item{\code{chd}}{a numeric vector event indicator for CHD at study
      exit (1=Yes, 0=No)}
    \item{\code{y}}{a numeric vector number of years at risk}
    \item{\code{hieng}}{a factor with levels \code{low} for energy <
      2750 or \code{high} for energy >= 2750}
    \item{\code{energy}}{a numeric vector total recorded energy intake (kCal/day/100)}
    \item{\code{job}}{a factor for occupation with levels for bus
      \code{driver}, bus \code{conductor} and \code{bank} work}
    \item{\code{month}}{a numeric vector month of the year (1-12)}
    \item{\code{height}}{a numeric vector height in cm}
    \item{\code{weight}}{a numeric vector weight in kg}
    \item{\code{doe}}{a Date for date of study entry}
    \item{\code{dox}}{a Date for date of study exit}
    \item{\code{dob}}{a Date for date of birth}
    \item{\code{yob}}{a numeric vector for continuous year of birth}
    \item{\code{yoe}}{a numeric vector for continuous year of entry}
    \item{\code{yox}}{a numeric vector for continuous year of exit}
  }
}
%% \details{
%% %%  ~~ If necessary, more details than the __description__ above ~~
%% }
%% \source{
%% %%  ~~ reference to a publication or URL from which the data were obtained ~~
%% }
%% \references{
%% %%  ~~ possibly secondary sources and usages ~~
%% }
\examples{
data(diet)
## maybe str(diet) ; plot(diet) ...
}
\keyword{datasets}
