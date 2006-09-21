\name{SPIKE_IN95}
\alias{SPIKE_IN95}
\docType{data}
\title{set of spike-in genes contained in the HGU95 dataset}
\description{This dataset contains a set of gene names contained in the HGU95 dataset}
\usage{data(SPIKE_IN95)}
\format{The set of spike-in gene names contained in the HGU dataset}
\source{This data is experimental data extracted from the publicly available HGU dataset}
\references{
 E. Carlon and T. Heim, Physica A 362, 433 (2006).
}
\seealso{
\code{\link{Ivsc}}, \code{\link{IvsDG}}, \code{\link{collapse}}, \code{\link{SPIKE_IN}}, \code{\link{hgu}}
, \code{\link{SPIKE_INA}}, \code{\link{SPIKE_INB}}, \code{\link{SPIKE_INH}}
}
\examples{
## you can first check if the data matches the predicted hybridisation value according to the langmuir 
## value, from the intensity versus the concentration value
Ivsc(SPIKE_IN95[1])

## you can then plot the value of the Intensity of the probe with the predicted value of the hybridisation
## according to the Delta G, value
IvsDG(SPIKE_IN95[4],128)

## The collapse function will finally plot all the values of the probe set according to 
## the langmuir absorption theory

collapse(SPIKE_IN95[2])

## By comparing the matched value and the mismatches, you will be able to identify errors which 
## could have done while sampling the data, or if the error happens repeatedly this will show errors 
## which will have happened while sequencing old data.
}
\keyword{datasets}
