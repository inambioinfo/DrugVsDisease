\name{profiles}
\alias{profiles}
\docType{data}
\title{
Gene Expression Profiles
}
\description{
List containing ranked lists of gene expression and associated p-values for a set of profiles.
}
\usage{data(profiles)}
\format{
	List containing rank of differential gene expression and pvalues in a list. Each item in the list contains a matrix. Matrix has number of rows equal number genes and columns for the different profiles.
}

\source{
\href{http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-22528}{http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-22528}
}

\examples{
data(profiles)
## maybe str(profiles) ; plot(profiles) ...
}
\keyword{datasets}
