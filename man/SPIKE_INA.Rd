\name{SPIKE_INA}
\alias{SPIKE_INA}
\docType{data}
\title{Artificial Spike-in probesets}
\description{This dataset contains the names of the probesets contained in the hgu dataset}
\usage{data(SPIKE_INA)}
\format{This dataset contains a set of String containing the names of the Artificial genes contained in the HGU dataset}
\references{
 E. Carlon and T. Heim, Physica A 362, 433 (2006).
}
\seealso{
\code{\link{Ivsc}}, \code{\link{IvsDG}}, \code{\link{collapse}}, \code{\link{SPIKE_IN}}, \code{\link{hgu}}
, \code{\link{SPIKE_INB}}, \code{\link{SPIKE_INH}}
}
\examples{
## you can first check if the data matches the predicted hybridisation value according to the langmuir 
## value, from the intensity versus the concentration value
Ivsc(SPIKE_INA[1])

## you can then plot the value of the Intensity of the probe with the predicted value of the hybridisation
## according to the Delta G, value
IvsDG(SPIKE_INA[4],128)

## The collapse function will finally plot all the values of the probe set according to 
## the langmuir absorption theory

collapse(SPIKE_INA[2])

## By comparing the matched value and the mismatches, you will be able to identify errors which 
## could have done while sampling the data, or if the error happens repeatedly this will show errors 
## which will have happened while sequencing old data.
}
\keyword{datasets}
