\name{Ivsc}
\alias{Ivsc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot of intensity vs. concentration for given probes }
\description{
     The function Ivsc plots intensity as a function of a concentration
     for a given probe in the spike-in Latin square experiments.  It also
     peforms a non-linear data fit (using the package nls in the R-package
     stats) of the experimental data using the Langmuir Isotherm:
              I = I0 + Ac/(K+c)
     Solid and dashed lines are best fits according to this formula.
     Imax in the plot are given by Imax=I0+A, ie the asymptotic
     intensity in the limit of c to infinity. 
}
\usage{
Ivsc(probe_set, probe = "NULL", outfile = "NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probe_set}{ Probe set number of the probe set analyzed }
  \item{probe}{ Integer giving the probe number (if not give the probe 1 is selected) }
  \item{outfile}{output the plotted data to a postscript file}
}
\references{
 E. Carlon and T. Heim, Physica A 362, 433 (2006).
}

\author{
  Delphine Baillon, Paul Leclercq, Sarah Ternisien, Thomas Heim and Enrico Carlon

Maintainer: Enrico Carlon <enrico.carlon@iemn.univ-lille1.fr>
}
\note{ eliminate output!
      If I-I0<0 set a small positive value instead (as 0.00001)
\section{Warning }{ Some probes have an irregular behavior and the non-linear square fit
     does not converge.}
}
\seealso{
  \code{\link{collapse}}, \code{\link{IvsDG}}, \code{\link{hgu}}, \code{\link{SPIKE_IN}}, \code{\link{SPIKE_INA}}, \code{\link{SPIKE_INB}}, \code{\link{SPIKE_INH}}
}
\examples{
     Ivsc("37777_at",4)
}
\keyword{ models }
\keyword{ data }
