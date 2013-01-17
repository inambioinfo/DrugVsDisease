classifyprofile <-
function(data,pvalues=NULL,case=c("disease","drug"),type=c("fixed","dynamic","range"),lengthtest=100,ranges=seq(100,2000,by=100),adj=c("qvalue","BH"),dynamic.fdr=0.05,signif.fdr=0.05,customRefDB=NULL, noperm=1000,customClusters=NULL,clustermethod=c("single","average"),avgstat=c("mean","median"),cytoout=FALSE,customsif=NULL,customedge=NULL,cytofile=NULL,no.signif=10){
	clustermethod<-match.arg(clustermethod)
	avgstat<-match.arg(avgstat)
	case<-match.arg(case)
	type<-match.arg(type)
	
	adj<-match.arg(adj)
	nodename=colnames(data)

	library(qvalue)
	if(!is.matrix(data)){
		data<-read.table(data,header=TRUE,row.names=1)
		data<-as.matrix(data)
	}
	if(!is.null(pvalues)&&!is.matrix(pvalues)){
		pvalues<-read.table(pvalues,header=TRUE,row.names=1)
		pvalues<-as.matrix(pvalues)
	}
#data<-as.matrix(data)
#pvalues<-as.matrix(pvalues)
	nodename=colnames(data)
	if(cytoout&&is.null(cytofile))stop('filename for cytoscape output needed')
	
	if(type=="dynamic"&&is.null(pvalues))stop('A numeric matrix of p-values or an lmobject is required to use p-values to select gene signatures')
	
	#read in customRefDB and customClusters if provided
	if(!is.null(customRefDB)&&!is.null(customClusters)){
		if(!is.matrix(customRefDB))customRefDB<-read.table(customRefDB,header=TRUE,row.names=1)
		if(!is.data.frame(customClusters))customClusters<-read.table(customClusters,header=TRUE)
		refnames<-colnames(customRefDB)
		refnames<-make.names(refnames)
		clustnames<-customClusters[,which(colnames(customClusters)%in%c("Drug","Disease"))]
		clustnames<-make.names(clustnames)
		intersection<-which(refnames%in%clustnames)
		if(length(refnames)!=length(intersection))stop('Profiles in custom rank profiles and clustom clusters do not match')
	}
		
	#calculate enrichment scores
	ESvals<-.calculateES(data,case=case,type=type,pvalues=pvalues,adj=adj,dynamic.fdr=dynamic.fdr,signif.fdr=signif.fdr,customRefDB=customRefDB,noperm=noperm,ranges=ranges,lengthtest=lengthtest)
	
			
	#assign significant enrichment scores to clusters		
	if(clustermethod=="single"){
			
		clusterassignments<-.classifysinglelinkage(ESvals$significance,customClusters=customClusters,case=case,no.signif=no.signif)
	}else{
		clusterassignments<-.classifyaveragelinkage(ESvals$scores,customClusters=customClusters,case=case,statistic=avgstat)
		
	}
	#check have at least one significant results
	nacheck<-unlist(clusterassignments)
	if(is.list(nacheck)){
		nacheck2<-unlist(nacheck)
		nas<-which(is.na(nacheck2))
		if(length(nas)==length(nacheck2)){stop('No significant matches found. Consider changing the FDR threshold')}

	}else{
		nas<-which(is.na(nacheck))
		if(length(nas)==length(nacheck)){stop('No significant matches found. Consider changing the FDR threshold')}
	}
	#generate cytoscape output
	if(cytoout){
		if(!is.null(customClusters)){
			if(is.null(customsif)|is.null(customedge)){
				stop('Custom SIF and Edge attribute file for cytoscape required since using custom clusters')
			}else{
				if(!is.data.frame(customsif)){customsif<-read.table(customsif,header=FALSE)}
				if(!is.data.frame(customedge)){customedge<-read.table(customedge,header=FALSE)}			
			}
		}
		.writecytoscape(case=case,clusterassignments=clusterassignments,filename=cytofile,customsif=customsif,customedge=customedge,nodename=nodename,customClusters=customClusters)
	}
	return(clusterassignments)
		
}
