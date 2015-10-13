\name{customClust}
\alias{customClust}
\docType{data}
\title{
Custom Clusters
}
\description{
Matrix containing a set of drugs and their corresponding clusters which can be used as input to the classifyprofile function.
}
\usage{data(customClust)}
\format{
Data frame with 47 rows and 2 columns. Columns are headed "Drug" and "Cluster". 
}
\details{

Each row refers to a compound in the customDB data set
\code{link{customDB}} with its corresponding cluster assignment in the
second column. These profiles are a subset of the Connectivity Map data
[1] (full set available in the DrugVsDiseasedata package
\pkg{DrugVsDiseasedata}, data object \code{drugRL}, for example use. Clusters
were generated using affinity propagation clustering [2]

}
\source{
\href{http://www.broadinstitute.org/cmap/}{http://www.broadinstitute.org/cmap/}
}
\references{
[1] Lamb J et~al. (2006) The Connectivity Map: Using Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease. Science, 313(5795), 1929-1935.
[2] U. Bodenhofer, A. Kothmeier, and S. Hochreiter. APCluster: an R package for affinity propagation clustering. Bioinformatics, 27(17):2463-2464, 2011.
}
\examples{
data(customClust)

}
\keyword{datasets}
