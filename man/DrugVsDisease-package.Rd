\name{DrugVsDisease-package}
\alias{DrugVsDisease-package}
\alias{DrugVsDisease}
\docType{package}
\title{
DrugVsDisease Package Overview
}
\description{
This package generates ranked lists of differential gene expression for either disease or drug profiles. Input data can be downloaded from Array Express [1] or GEO [2], or from local CEL files. Ranked lists of differential expression and associated p-values are calculated using Limma [3]. Enrichment scores [4] are calculated to a reference set of default drug or disease profiles, or a set of custom data supplied by the user. Significance scores are output in Cytoscape \href{http://www.cytoscape.org/}{http://www.cytoscape.org/} format.
}
\details{
\tabular{ll}{
Package: \tab DvD\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-06-15\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
Profiles are calculated via generateprofiles and selected profiles are then classified using classifyprofile.
}
\author{
C. Pacini

Maintainer: J Saez-Rodriguez <saezrodriguez@ebi.ac.uk>
}
\references{
[1]Parkinson et al. (2010) ArrayExpress update an archive of microarray and high-throughput sequencing-based functional genomics experiments. Nucl. Acids Res.,doi: 10.1093/nar/gkq1040.

[2]Barrett T et al. (2011) NCBI GEO: archive for functional genomics data sets-10 years on. Nucl. Acids Res, 39, D1005-D1010.

[3]Smyth et al. (2004). Linear models and empirical Bayes method for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology, Vol. 3, No. 1, Article 3. 

[4]Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L.,
Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander,
E. S. & Mesirov, J. P. (2005) \emph{ Gene set enrichment analysis: A
knowledge-based approach for interpreting genome-wide expression profiles.}
Proc. Natl. Acad. Sci. USA 102, 15545-15550.


}

\keyword{ package }
\seealso{
\code{\link{generateprofiles}}, \code{\link{selectrankedlists}}, \code{\link{classifyprofile}}
}
\examples{
profiles<-generateprofiles(input="GEO",accession="GDS2617",case="disease",statistic="t",annotation="hgu133a")
selprofile<-selectrankedlists(profiles,6)
classification<-classifyprofile(data=selprofile$ranklist,signif.fdr=1,noperm=20)}

