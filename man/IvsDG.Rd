\name{IvsDG}
\alias{IvsDG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot Intensity as function of the affinity for a given probe set at
fixed concentration. }
\description{
     IvsDG plots intensity vs affinity (or free energy) for a probe
     set at a given concentration. The outcome is compared with the
     prediction from the Langmuir isotherm at that concentration.
     Two graphs are shown: on the left intensity vs. probe number for PM
     (blue) and MM (red); on the right the same value plotted as function
     of the affinities. The black line is the Langmuir Isotherm at the
     given concentration. The two green lines correspond to concentrations
     fourfold higher and lower compared to the given one.
}
\usage{
IvsDG(probe_set, conc, outfile = "NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probe_set}{ Probe set number of the probe set analyzed }
  \item{conc}{ Concentration value }
  \item{outfile}{ "PS" output on a postscript file }
}
\references{
 E. Carlon and T. Heim, Physica A 362, 433 (2006).
}

\author{
  Delphine Baillon, Paul Leclercq, Sarah Ternisien, Thomas Heim and Enrico Carlon

Maintainer: Enrico Carlon <enrico.carlon@iemn.univ-lille1.fr>
}
\seealso{
  \code{\link{Ivsc}}, \code{\link{collapse}}, \code{\link{hgu}}, \code{\link{SPIKE_IN}}, \code{\link{SPIKE_INA}}, \code{\link{SPIKE_INB}}, \code{\link{SPIKE_INH}}
}
\examples{
     data(hgu)
     IvsDG("1024_at",64)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ data }% __ONLY ONE__ keyword per line
