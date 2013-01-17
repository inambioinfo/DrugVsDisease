\name{generateprofiles}
\alias{generateprofiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Profiles
}
\description{
Processing Affymetrix data to generate ranked lists of differential gene expression and associated p-values. 
}
\usage{
generateprofiles(input = c("AE", "GEO", "localAE", "local"), normalisation = c("rma", "mas5"), accession = NULL, customfile = NULL, celfilepath = NULL, sdrfpath = NULL, case = c("disease", "drug"), statistic = c("coef", "t", "diff"), annotation = NULL, factorvalue = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{input}{
Character string denoting the source of the data. One of AE (default), GEO, localAE or local.
}
  \item{normalisation}{
Character string denoting the normalisation procedure as implemented in the affy package. One of mas5 (default) or rma.
}
  \item{accession}{
Optional character string giving the database reference for use with either the AE or GEO options.
}
  \item{customfile}{
Optional character string giving the path of a file containing the factor values associated with the CEL files specified in folder celfilepath
}
  \item{celfilepath}{
Optional character string giving the path of a folder containing CEL files to analyse.
}
  \item{sdrfpath}{
Optional character string giving path of an sdrf file corresponding to CEL files in celfilepath
}
  \item{case}{
Character string, one of disease (default) or drug denoting whether the input profiles are disease or drug profiles.
}
  \item{statistic}{
Character string, one of coef (default), t or diff. 
}
  \item{annotation}{
Optional character string giving the platform of the affymetrix files
}
  \item{factorvalue}{
Optional character string giving the name of the factor value in the GEO database.
}
}
\details{
Input types of AE and GEO use raw data download from Array Express using the ArrayExpress [1] package or processed GDS files from GEO using the GEOquery package [2]. CEL files and sdrf files downloaded from Array Express and stored locally can be processed using localAE option with the sdrf file path specified in sdrfpath and the path of the folder containing the CEL files contained in celfilepath. Users data stored locally can be processed using the local option with CEL file folders in celfilepath and factors associated with the CEL files in customfile. Where metadata may be missing from the GEO database, platform annotations can be specified using the annotation parameters and the name of main factor value (e.g. disease status, or compound treatment) using factorvalue option. 
Raw CEL files are normalised (rma or mast)[3] and data is converted from probes to genes using BioMart annotations [4]. Linear models are fitted using the database factor vales or user provided factors for locally stored data [5]. The differential expression is calculated for HUGO genes with the mapping performed automatically for Affymetrix platforms, HGU133A, HGU133Plus2 and HGU133A2 using BioMart. The differential expression statistic is one of coef (default), which corresponds to log (base 2) FC, diff (which is the difference between raw (non-logged) expression values, or t for the t-statistic based on log base 2 expression values.
}
\value{
List with two elements:
\item{Ranklist}{Matrix containing the ranks of gene expression. Rows containing the genes, columns the different profiles}
\item{Pvalues}{Matrix containing the associated p-values to the differential expression profiles in Ranklist}
}


\references{
[1]Kauffmann et al. (2009) Importing Array Express datasets into R/Bioconductor. Bioinformatics, 25(16):2092-4.

[2]Davis et al. (2007) GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. Bioinformatics, 14, 1846-1847.

[3]Irizarry et al. (2003) Summaries of Affymetrix GeneChip probe level data. Nucleic Acids Research, 31(4); e15.


[4]Durinck et al. (2009). Mapping identifiers for the integration of genomic datasets with the R/Bioconductor package biomaRt. Nature Protocols  4, 1184-1191.

[5]Smyth et al. (2004). Linear models and empirical Bayes method for assessing differential expression in microarray experiments. Statistical Applications in Genetics and Molecular Biology, Vol. 3, No. 1, Article 3. 
}
\author{
C. Pacini
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{classifyprofile}}
}
\examples{

profiles<-generateprofiles(input="GEO",accession="GDS2617",case="disease",statistic="t",annotation="hgu133a")


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
