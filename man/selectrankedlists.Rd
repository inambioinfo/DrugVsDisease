\name{selectrankedlists}
\alias{selectrankedlists}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select Ranked Lists
}
\description{
Given a list with two elements, one containing ranklists for a set of regression models and the second containing the associated p-values. This function is used to extract a subset of the models.
}
\usage{
selectrankedlists(ranklist, colsinc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ranklist}{
List with two elements (as output from generate profiles):
\emph{Ranklist}: Matrix containing the ranks of gene expression. Rows containing the genes, columns the different profiles.

\emph{Pvalues}: Matrix containing the associated p-values to the differential expression profiles in Ranklist.
}
 \item{colsinc}{
Vector of integers containing the column references of the profiles to select.
}
}
\details{
Format of list provided to selectrankedlists is the same as is output from generateprofiles \code{\link{generateprofiles}}. The output from selectrankedlists can be used as input to the classify profile function \code{\link{classifyprofile}}.
}
\value{
\item{Ranklist}{Matrix containing the ranks of gene expression. Rows containing the genes, columns the selected profiles}
\item{Pvalues}{Matrix containing the associated p-values to the differential expression profiles in Ranklist}

}

\author{
C. Pacini
}

\source{
\href{http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-22528}{http://www.ebi.ac.uk/arrayexpress/experiments/E-GEOD-22528}
}
\seealso{
\code{\link{generateprofiles}}, \code{\link{classifyprofile}}.
}
\examples{
data(profiles)
selectprofile<-selectrankedlists(profiles,1)
classification<-classifyprofile(data=selectprofile$ranklist,noperm=10,signif.fdr=1)}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
