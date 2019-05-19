\name{Simulating Clones}
\alias{generateTumorData}
\alias{plotTumorData}
\alias{tumorGen}
\alias{dataGen}
\title{ Simulating Tumor Clones }
\description{
  Simulating copy number segmentation data and sequencing mutation data
  for tumors composed of multiple clones.
}
\usage{
generateTumorData(tumor, snps.seq, snps.cgh, mu, sigma.reads,
                  sigma0.lrr, sigma0.baf, density.sigma)
plotTumorData(tumor, data)
tumorGen(...)
dataGen(tumor, ...)
}
\arguments{
  \item{tumor}{ an object of the \code{Tumor} class. }
  \item{snps.seq}{an integer; the total number of germline variants and
  somatic mutations to simulate in the tumor genome.}
  \item{snps.cgh}{an integer; the number of single nucleotide
    polymorphisms (SNPs) to simulate as meaurements made to estimate
    copy number.}
  \item{mu}{an integer; the average read depth of a simulated sequencing
    study giving rise to mutations.}
  \item{sigma.reads}{a real number; the standard deviation of the number
    of simulated sequencing reads per base.}
  \item{sigma0.lrr}{a real number; the standard deviation of the
    simulated per-SNP log R ratio (LRR) for assessing copy number.}
  \item{sigma0.baf}{a real number; the standard deviation of the
    simulated B allele frequency (BAF) for assessing copy number.}
  \item{density.sigma}{a real number; the standard deviation of a beta
    distribution used to simulate the number of SNP markers per copy
    number segment.}
  \item{data}{ a list containing two data frames, \code{cn.data} and
    \code{seq.data}, as produced by \code{generateTumorData}.}
  \item{...}{additional variables}
}
\details{
  Copy number and mutation data are simulated essentially
  independently. Each simulation starts with a single "normal" genome, 
  and CNVs and/or mutations are randomly generated for each new "branch"
  or subclone. (The number of subclones depends on the input
  parameters.) Each successive branch is randomly determined to descend
  from one of the existing clones, and therefore contains both the
  aberrations belonging to its parent clone and the novel aberrations
  assigned to it. Depending on input parameters, the algorithm can also
  randomly select some clones for extinction in the process of
  generating the heterogeneous tumor, to yield a more realistic
  population structure. 

  Note that \code{tumorGen} (an alias for \code{Tumor} that returns a
  list instead of a Tumor object) and \code{dataGen} (an alias for
  \code{generateTumorData}) are DEPRECATED.
}
\value{
  The \code{generateTumorData} function returns a list with two
  components, \code{cn.data} and \code{seq.data}. Each component is
  itself a data frame. Note that in some cases, one of these data frames
  may have zero rows or may be returned as an \code{NA}.

  The \code{cn.data} component contains seven columns:
  \describe{
    \item{\code{chr}}{the chromosome number;}
    \item{\code{seq}}{a unique segment identifier;}
    \item{\code{LRR}}{simulated segment-wise log ratios;}
    \item{\code{BAF}}{simulated segment-wise B allele frequencies;}
    \item{\code{X} and \code{Y}}{simulated intensities for two separate
      alleles/haplotypes per segment; and} 
    \item{\code{markers}}{the simulated number of SNPS per segment.}
  }

  The \code{seq.data} component contains eight columns:
  \describe{
    \item{\code{chr}}{the chromosome number;}
    \item{\code{seq}}{a unique "segment" identifier;}
    \item{\code{mut.id}}{a unique mutation identifier;}
    \item{\code{refCounts} and \code{varCounts}}{the simulated numbers
      of reference and variant counts per mutation;}
    \item{\code{VAF}}{the simulated variant allele frequency;}
    \item{\code{totalCounts}}{the simulated total number of read counts;
      and}
    \item{\code{status}}{a character (that should probably be a factor)
      indicating whether a variant should be viewed as somatic or
      germline.}
  }

  The \code{plotTumorData} function invisibly returns its \code{data}
  argument.
}
\author{
  Kevin R. Coombes \email{krc@silicovore.com},
  Mark Zucker \email{zucker.64@buckeyemail.osu.edu}
}
\examples{
psis <- c(0.6, 0.3, 0.1) # three clones
# create tumor with copy number but no mutation data
tumor <- Tumor(psis, rounds = 400, nu = 0, pcnv = 1, norm.contam = FALSE)
# simulate the dataset
dataset <- generateTumorData(tumor, 10000, 600000, 70, 25, 0.15, 0.03, 0.1)
#plot it
plotTumorData(tumor, dataset)
}
\keyword{array}
