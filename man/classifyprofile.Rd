\name{classifyprofile}
\alias{classifyprofile}
\title{
Classify Profiles
}
\description{
For a set of ranked gene expression profiles, enrichment scores to a default or custom set of profiles. The ranked gene expression profiles can have been generated from generate profiles or provided by the user. The enrichment scores are assessed for significance using permutations to generate random profiles.
}
\usage{
classifyprofile(data, pvalues = NULL, case = c("disease", "drug"), type = c("fixed", "dynamic", "range"), lengthtest = 100, ranges = seq(100, 2000, by = 100), adj = c("BH","qvalue"), dynamic.fdr = 0.05, signif.fdr = 0.05, customRefDB = NULL, noperm = 1000, customClusters = NULL, clustermethod = c("single", "average"), avgstat = c("mean", "median"), cytoout = FALSE, customsif = NULL, customedge = NULL, cytofile = NULL, no.signif = 10,stat=c("KS","WSR"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix gene expression profiles. Rows are genes whose names match the genelist from DvDdata and columns are the different profiles. This can also be a path to a txt file containing the data.
}
  \item{pvalues}{
Optional numerical matrix of pvalues associated with the differential expression of the ranked lists in data argument. This can also be a path to a txt file containing the data.
}
  \item{case}{
Character string indicating whether or not the input profiles are disease (default) or drug profiles
}
  \item{type}{
Character string giving the method of gene set sizes, one of fixed (default), dynamic or range
}
  \item{lengthtest}{
Integer giving the number of genes to generate a set size for use with type=fixed, default 100.
}
  \item{ranges}{
A vector of integer values giving a range of gene set sizes, for use with type=range, default 100 to 2000 every 100.
}
  \item{adj}{
Character string to set the method for multiple hypothesis testing, one of qvalue or BH (default)
}
  \item{dynamic.fdr}{
Double giving the false discovery rate to determine the gene set sizes (default 5%) 
}
  \item{signif.fdr}{
Double giving the false discovery rate to determine significance of enrichment scores (default 5%)
}
  \item{customRefDB}{
Optional matrix of reference ranked profiles to compare the input profiles to. Alternatively a string giving the name of the path where the database is stored.
}
  \item{noperm}{
Integer for the number of permutation profiles (default 1000)
}
  \item{customClusters}{
Optional data frame of cluster assignments which relate to customRefDB
}
  \item{clustermethod}{
Character string to give the cluster method one of single (default) or average
}
  \item{avgstat}{
Character string for the statistic used with average cluster method. One of mean (default) or median.
}
  \item{cytoout}{
Logical if Cytoscape SIF and Edge Attribute files should be produced (default is FALSE)
}
  \item{customsif}{
Optional SIF input needed for use with custom clusters. Can be an R object or a character string containing a file path
}
  \item{customedge}{
Optional Edge attribute input needed for use with custom clusters. Can be an R object or a character string containing a file path
}
  \item{cytofile}{
Character string for the filename of the Cytoscape output
}
  \item{no.signif}{
Integer giving the maximum number of significant enrichment scores to return. Default is 10.
}
\item{stat}{One of KS (default) or WSR. KS is the Kolmogorov-Smirnov type statistic which equally weights all elements of the gene set. WSR uses both the sign and the position in the ranked list to calculate an enrichment score.}
}
\details{
The classify profile function contains a default set of drug (from the Connectivity Map [2]) and disease profiles (from various GEO profiles) with corresponding clusters which input profiles of disease and drug respectively are compared to. Enrichment scores [1] are calculated between the input profiles and the corresponding inverted reference profiles such that, the score measures the enrichment of up-regulated genes in the input profiles in the down-regulated genes in the reference set[3]. The gene set sizes to use in the Enrichment scores can be specified using one of three methods - fixed, range or dynamic. With the fixed method the user specifies a fixed number of genes to use with all input profiles. The range option takes a vector of integers - enrichment scores are calculated for each profile using a gene set size as given by the vector. The dynamic option uses the p-values to determine the gene set size according to the number of significantly differential expressed genes (following multiple hypothesis correction). Multiple hypothesis correction is done using one of two methods, qvalue or Benjamini-Hochberg.
Two enrichment scores are calculated for the up and down regulated scores. These contain the scores where the input profiles gene set is compared to the reference profile and where the reference profiles equivalent gene set (as determined by the type method) are compared to the input profile when using the KS option [4]. The WSR option implements the score of [5]. Various optional parameters exist for comparing an input profile to a users own reference set. For using custom reference data the user needs to provide the custom ranked lists of differential expression (customRefDB), the corresponding clusters (or network) between nodes in the reference data set (customClusters) and the SIF and Edge attribute files for Cytoscape option if cytoout=TRUE. Default clusters are provided by the DvDdata package, and include a drug and disease network. Input profiles to classifyprofile can be assigned to clusters using either single or average linkage. With single linkage an edge is drawn between the input profile and any significant scoring (up to a user defined maximum of no.signif) reference profiles. For average linkage either the mean or median (specified through avgstat) of the scores to each member in a cluster is calculated and the profile is assigned to the cluster with the highest average score.

To use your own preprocessed data, make sure the txt files for the data (and optional pvalues) have rownames with genes matching those in the reference data set. The files should have genes as row names in the first column and the header (col names) giving the names of the input profile(s). The input to classifyprofile is then a string of the path to the files.
}
\value{
Data Frame for each profile with elements:
\item{Node}{Names of the profiles with significant scores to the input profile(s)}
\item{ES Distance}{The distance between the input profile and corresponding reference profile. Defined as 1-ES}
\item{Cluster}{The cluster number of the node}
\item{RPS}{Running sum Peak Sign, taking values -1 to indicate an inverse relationship (potentially) therapeutic, or 1 to indicate a similar profile.}
}
\references{
[1]Subramanian A et~al. (2005) Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. PNAS, 102(43), 15545-15550.

[2]Lamb J et~al. (2006) The Connectivity Map: Using Gene-Expression Signatures to Connect Small Molecules, Genes, and Disease. Science, 313(5795), 1929-1935.

[3]Sirota M et~al. (2011) Discovery and Preclinical Validation of Drug Indications Using Compendia of Public Gene Expression Data. Sci Trans Med,3:96ra77.

[4]Iorio et al. (2010) Discovery of drug mode of action and drug repositioning from transcriptional responses. PNAS, 107(33), 14621-14626.

[5]Zhang S et al. (2008) A simple and robust method for connecting small- molecule drugs using gene-expression signatures. BMC Bioinformatics, 9:258.
}
\author{
C. Pacini
}




\seealso{
Function for generating profiles for input to classifyprofile: \code{\link{generateprofiles}}.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(selprofile)
classification<-classifyprofile(data=selprofile$ranklist,signif.fdr=1,noperm=20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{classify}
\keyword{profile}
\keyword{GSEA}
