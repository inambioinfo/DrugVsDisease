\name{combineProfiles}
\alias{combineProfiles}
\title{
Combine Profiles
}
\description{
Combine a set of hybridisations using the median rank method as proposed by Warnat
}
\usage{
combineProfiles(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A list of R objects, each object should contain a matrix of expression values to be merged with each other. The first object will be used as the reference data set.
}
}
\details{
The combine profiles function merges a set of expression values, normalising across different experiments to facilitate the meta-analysis of similar experiments i.e. cell lines treated with the same drug but which may be from different "batches" or platforms.}
\value{
Matrix of rank normalised expression data. Rows are genes, columns are hybridisations from the different experiments.}
\references{
R code from CONOR package. 
Need also reference to the method.
}
\author{
C. Pacini
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classify}
\keyword{profile}
\keyword{GSEA}
