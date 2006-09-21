\name{SPIKE_INH}
\alias{SPIKE_INH}
\docType{data}
\title{Human Spike-in probe-set names}
\description{
   This dataset contains the names of the Human probe-sets contained in the HGU dataset
}
\usage{data(SPIKE_INH)}
\format{names of the human probe-sets contained in the HGU dataset}
\references{
 E. Carlon and T. Heim, Physica A 362, 433 (2006).
}
\seealso{
\code{\link{Ivsc}}, \code{\link{IvsDG}}, \code{\link{collapse}}, \code{\link{SPIKE_IN}}, \code{\link{hgu}}
, \code{\link{SPIKE_INA}}, \code{\link{SPIKE_INB}}
}
\examples{
## you can first check if the data matches the predicted hybridisation value according to the langmuir 
## value, from the intensity versus the concentration value
Ivsc(SPIKE_INH[3])

## you can then plot the value of the Intensity of the probe with the predicted value of the hybridisation
## according to the Delta G, value
IvsDG(SPIKE_INH[5],256)

## The collapse function will finally plot all the values of the probe set according to 
## the langmuir absorption theory

collapse(SPIKE_INH[2])

## By comparing the matched value and the mismatches, you will be able to identify errors which 
## could have done while sampling the data, or if the error happens repeatedly this will show errors 
## which will have happened while sequencing old data.
}
\keyword{datasets}
