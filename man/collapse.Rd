\name{collapse}
\alias{collapse}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Data collapse of all concentrations into a single graph }
\description{    This function takes as input one or more (up to four) probe sets
     of the Latin square spike-in data and produces collapse plots.
     A collapse plot contains data of different concentrations into
     a single graph. The user can compare in how far the data follows
     the predicted Langmuir behavior which is also given in the plot.
     Two models are compared: the basic Langmuir Isotherm and the Langmuir
     Isotherm with hybridization in solution.
}
\usage{
collapse(probe_set, param = "NULL", probes = "NULL", output = "NULL", filename = "NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probe_set}{ This has to take the value of a probe set }
  \item{param}{ In input one or more probe sets can be given }
  \item{probes}{ A vector containing the probes  }
  \item{output}{ "PS" output on a postscript file }
  \item{filename}{ the file in which collapses are given }
}
\references{ E. Carlon and T. Heim, Physica A 362, 433 (2006).}
\seealso{
  \code{\link{Ivsc}}, \code{\link{IvsDG}}, \code{\link{hgu}}, \code{\link{SPIKE_IN}}, \code{\link{SPIKE_INA}}, \code{\link{SPIKE_INB}}, \code{\link{SPIKE_INH}}
}
\author{   Delphine Baillon, Paul Leclercq, Sarah Ternisien, Thomas Heim and Enrico Carlon }
\examples{
## You may display the matched intensities of a Probe-Set according to the Delta-G value
     collapse("1091_at")

## You may restrict the value to the Perfect match or mis-matches
     collapse("1091_at","PM")

## You may restrict the values risplayed for only a number of probes
     collapse("1091_at",probes=c(1,9))

## You may output the graphs to a postscript file
     collapse("1091_at",output="PS",filename="outfile.ps")

## You may display up to 4 probe-sets in the same window
     collapse(c("1091_at","37777_at",SPIKE_INA[1:2]))

## You can also use the values of the probe-sets contained in one of the Vectors of Human, Bacteria, 
## or Artificial Probe-sets
     collapse(SPIKE_INH)
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ data }% __ONLY ONE__ keyword per line
