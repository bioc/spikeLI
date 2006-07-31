\name{hgu}
\alias{hgu}
\docType{data}
\title{Selected Probe Set data}
\description{
 This selected probe sets information contains the sequence of the selected probe sets, as well as the match and mismatch onformation and Delta G value required for the langmuir analysis
}
\usage{data(hgu)}
\format{
  A data frame with 11452 observations on the following 9 variables.
  \describe{
    \item{\code{Probe.Set.Name}}{Name of probe set}
    \item{\code{conc}}{a numeric vector}
    \item{\code{Ipm}}{a numeric vector}
    \item{\code{Imm}}{a numeric vector}
    \item{\code{Seq}}{DNA Sequence of the probe}
    \item{\code{DGpm}}{DG value of perfect match of the probe}
    \item{\code{DGmm}}{Delta G value of the mismatch of the probe}
    \item{\code{DGRNA}}{Delta G value of the RNA}
    \item{\code{FILE}}{a factor with levels \code{HGU133} \code{HGU95}}
  }
}
\references{
 E. Carlon and T. Heim, Physica A 362, 433 (2006).
}  
\seealso{
\code{\link{Ivsc}}, \code{\link{IvsDG}}, \code{\link{collapse}}, \code{\link{SPIKE_IN}}, \code{\link{SPIKE_INA}}, \code{\link{SPIKE_INB}}, \code{\link{SPIKE_INH}}
}
\examples{
data(hgu)
## maybe str(hgu) ; plot(hgu) ...
}
\keyword{datasets}
