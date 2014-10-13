\name{selprofile}
\alias{selprofile}
\docType{data}
\title{
List: Differential gene expression and p-values
}
\description{
Ranklist which has the rank of each gene according to differential expression. P-values is the second element in the list containing the associated p-value for gene differential expression.
}
\usage{data(selprofile)}
\format{
List with two elements, the first ranklist containing the ranked position of the gene according to its differential expression. The second item in the list is the associated pvalues.
}
\details{
Example ranked list of differential expression. Taken from Array express experiment E-GEOD-22528
}
\source{
\href{http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-22528}{http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-22528}
}

\examples{
data(selprofile)

## maybe str(selprofile) ; plot(selprofile) ...
}
\keyword{datasets}
