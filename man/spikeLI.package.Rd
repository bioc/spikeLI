\name{spikeLI-package}
\alias{spikeLI-package}
\alias{spikeLI}
\docType{package}
\title{
   Analysis of Affymetrix spike-in data (HGU95 and HGU133 Latin square) 
   using the Langmuir Isotherm.
}
\description{
   spikeLI performs a series of analysis of Affymetrix spike-in data
   using inputs from physical-chemistry. It illustrates the advantages
   of such approach in determining expression levels and in identifying
   outliers compared to other methods. The analysis so far is restricted
   to spike-in genes. It will be extended to a generic CEL file.
   spikeLI does not require affy (and it is independent of any other
   bioconductor packages) as it reads spike-in data from a data frame
   variable hgu which is contained in the package.
}
\details{
\tabular{ll}{
Package: \tab spikeLI\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-05-05\cr
License: \tab GNU Public License\cr
}
    The package contains three basic functions:
    - Ivsc plot intensities as function of spike-in concentration for a
    fixed probe.
    - IvsDG plot intensities as function of affinity for a given probe
    set at fixed concentration.
    - collapse plot of intensities both as a function of concentration
    and affinities.
}

\author{
  Delphine Baillon, Paul Leclercq, Sarah Ternisien, Thomas Heim and Enrico Carlon

Maintainer: Enrico Carlon <enrico.carlon@polytech-lille.fr>
}
\seealso{
  \code{\link{collapse}}, \code{\link{Ivsc}}, \code{\link{IvsDG}}, \code{\link{hgu}}, \code{\link{SPIKE_IN}}, \code{\link{SPIKE_INA}}, \code{\link{SPIKE_INB}}, \code{\link{SPIKE_INH}}
, \code{\link{SPIKE_IN95}}
}
\references{
 E. Carlon and T. Heim, Physica A 362, 433 (2006).
}
\keyword{ package }
\keyword{ models }
