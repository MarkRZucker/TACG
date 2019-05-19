\name{Tumor-class}
\alias{Tumor-class}
\alias{Tumor}
\alias{getClone}
\alias{summary,Tumor-method}
\alias{coerce,Tumor,list-method}
\docType{class}
\title{Class "Tumor"}
\description{
  A class that represents tumors, thought of as a collection of
  (sub)clones each with an associated measure as a fraction of all
  tumor cells.
}
\usage{
Tumor(psi, rounds, nu = 100, pcnv = 0.5, norm.contam = FALSE, cnmax = 4)
getClone(tumor, i)
}
\arguments{
  \item{psi}{ a numeric vector containing non-negative values.}
  \item{rounds}{an integer; the number of generations through which to
    evolve the potential clones.}
  \item{nu}{an integer; the expected number of mutations in each clonal
    generation.}
  \item{pcnv}{a real number between 0 and 1; the probability of a CNV
    occuring per generation.}
  \item{norm.contam}{a logical value; should we treat one of the cell
    populations as normal cells that are "contaminating" the tumor
    specimen?}
  \item{cnmax}{an integer, the maximum copy number allowed in the
    simulated data.}
  \item{tumor}{ an object of the \code{Tumor} class.}
  \item{i}{a integer; which clone to extract.}
}
\details{
  The \code{Tumor} class is used to represent compex tumors, each of
  which consists a set of subclones representing different fractional
  parts of the tumor. Each clone is characterized by a set of copy
  number variants (modeled by the output produced by something like the
  \code{DNAcopy} package) and, optionally, a set of sequence
  mutations. Each of these genetic events is mapped to a specific
  interval or point in the human genome.
  
  In the current implementation, a \code{Tumor} consists of a weight
  vector that specifies the fractions of cells for each clone and a list
  of clones. At present, each clone is itself a list containing one (if
  there are no mutations) or two (if there are both copy number
  variants and mutations) data frames. This structure is likely to
  change in later versions of the package, since we expect to implement a
  full-fledged S4 class to represent clones. So, one should not rely on
  the current implementation. 
}
\value{
  The constructor returns a valid object of the \code{Tumor}
  class.
}
\section{Objects from the Class}{
  Although objects can be created using \code{new}, the preferred method
  is to use the constructor function, \code{Tumor}.
}
\section{Slots}{
  \describe{
    \item{\code{psi}}{ a \code{\link{WeightVector}} containing
      non-negative values whose sum  equals one. }
    \item{\code{clones}}{ a list, each of whose elements represents a
      clone.}
  }
}
\section{Methods}{
  \describe{
    \item{coerce(from, to, strict = TRUE)}{Convert the
      \code{Tumor} object into a simple numeric vector. Never
      actually used in this form, since the preferred method is to write
      \code{as(WV, "list")}.}
  }
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Mark Zucker \email{zucker.64@buckeyemail.osu.edu}
}
\seealso{
  \code{\link{WeightVector}}
}
\examples{
showClass("Tumor")

tumor <- Tumor(c(0.5, 0.3, 0.2), 100)
simpleTumor <- as(tumor, "list")
reformed <- as(simpleTumor, "Tumor")
}
\keyword{classes}

