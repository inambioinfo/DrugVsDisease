.datafromGEO <-
function(accession,destdir=".",factorvalue=NULL,case=c("disease","drug")){
	case<-match.arg(case)

	#make sure accession starts with GDS…
	testname<-substr(accession,1,3)
	if(testname!="GDS"){stop('Accession must be of type GDS, check accession name')}
	#if have a GDS accession start reference, try to download file from GEO:
	if(file.exists(paste(destdir,"/",accession,".soft.gz",sep=""))){file.remove(paste(destdir,"/",accession,".soft.gz",sep=""))}
	data<-getGEO(accession,destdir=destdir)
	#convert preprocessed data into expression set.
	eset<-GDS2eSet(data,do.log2=TRUE)
	#extract the names of the factors available in the annotation of the GDS data
	factornames<-names(pData(eset))
	platform<-annotation(eset)
	#check that we have a factor value to use in regression
	if(is.null(factorvalue)){
		sel<-which(factornames%in%c("disease state","dose","disease.state"))
	}else{
		sel<-which(factornames==factorvalue)
	}
	
	if(length(sel)==0&&length(factornames)!=0){
		fn<-c()
		for(i in 1:length(factornames)){
			fn<-paste(fn,factornames[i])
			
		}
		
		stop(paste('Neither disease nor dose/drug variable available for data set. Available factors: ',fn,' use parameter factorvalue to select a factor'))
	}
	
	if(length(sel)==0&&length(factornames)==0){stop('No factornames available in the pData of the eset. Consider processing using the local option.')}
	
	#load the list of factor values used to annotate experiments in the GEO database
	data(GEOfactorvalues,package="DrugVsDiseasedata")
	#extract those factors which are in the list of GEO factors
	factors<-GEOfactorvalues[which(GEOfactorvalues%in%factornames)]
	customfactors<-data.frame(pData(eset)[factors])
	CEL<-rownames(customfactors)
	customfactors<-data.frame(customfactors,CEL)
	factors<-colnames(customfactors)
	#find the factor (default or selected by user with parameter factorvalue) to set as either the "disease" or the "compound" factor - this is the label expected by function fitlm. This factor is used as the main explanatory variable to regress on in the fitlm function. 
	if(case=="disease"){
		if(is.null(factorvalue)){
			sel<-which(factors%in%c("disease state","disease.state"))
			
		}else{
			sel<-which(factors==factorvalue)
	
		}
		colnames(customfactors)[sel]<-"disease"
	}else{
		if(is.null(factorvalue)){
			sel<-which(factors=="dose")
			
		}else{
			sel<-which(factors==factorvalue)
			
		}
		colnames(customfactors)[sel]<-"compound"
	}
	file.remove(paste(destdir,"/",accession,".soft.gz",sep=""))
	return(list(normalised=exprs(eset),customfactors=customfactors,platform=platform))
}
.readlocalCEL <-
function(path){
	#assume reading affymetrix CEL files from a local directory
	affyobject<-ReadAffy(celfile.path=path)
	return(affyobject)
}
.readlocalAE <-
function(celpath,sdrfpath){
	#read CEL files from AE and corresponding sdrf file
	try(CELfiles<-ReadAffy(celfile.path=celpath))
	if(class(CELfiles)[1]=="try-error"){stop(paste('Unable to read CEL files from path: ',celpath,'. Check CEL files and path directory.',sep=""))}
	sdrffile<-read.delim(sdrfpath)
	filenames<-sdrffile[,"Array.Data.File"]
	sortsdrf<-sdrffile[order(sdrffile[,"Array.Data.File"]),]
	#now have sdrf rows matching the order of the CEL files
 	return(list(affyobject=CELfiles,factors=sortsdrf))
}
.datafromAE <-
function(experiment){

	#Get the data from Array Express package and factors from the sdrf file
	try(experimentAEset<-ArrayExpress(experiment))
	if(class(experimentAEset)[1]=="try-error"){stop(paste("Unable to download experiment",experiment,"from Array Express. Please check accession reference and platform for files. Consider using options 'local' or 'localAE' instead."))}
	platform<-annotation(experimentAEset)
	return(list(affyobject=experimentAEset,factors=pData(experimentAEset),platform=platform))		
		
	
}
.normalisedata <-
function(eset,normalisation=c("rma","mas5")){
	normalisation<-match.arg(normalisation)
	if(normalisation=="rma"){
			normalised<-rma(eset)
		
	}else{
		try(normalised<-mas5(eset))
		if(class(normalised)[1]=="try-error"){stop('Error in mas5, check annotation of CEL files are the same')}
		
	}
return(normalised)
}
.calculateES <-
function(data,type=c("fixed","dynamic","range"),adj=c("qvalue","BH"),lengthtest=250,ranges=seq(100,2000,by=100),pvalues=NULL,dynamic.fdr=0.01,signif.fdr=0.01,customRefDB=NULL,case=c("disease","drug"),noperm=5000,stat=c("KS","WSR")){
	case=match.arg(case)
	type=match.arg(type)
	
  	adj=match.arg(adj)
	experimentnames<-colnames(as.matrix(data))
	#data checks here
	if(!is.null(data)&&!is.numeric(data)){stop('Data must be a numeric matrix of expression values or rank numbers')}
	
	#if a customRefDB is not supplied load the default profiles from the cMap2data or DrugVsDiseasedata packages depending on whether the input profile is a drug or disease profile.
	if(is.null(customRefDB)){
		if(case=="disease"){
			data(drugRL,package="cMap2data")
			#drugRL is ranks in decreasing order so rank 1 is highest increase DE
			refDB<-drugRL
			if(stat=="WSR"){
				nogenes<-nrow(refDB)
				splitno<-round(nogenes/2)+1
				refDB<-splitno-refDB
			}
		}else{
			data(diseaseRL,package="DrugVsDiseasedata")
			#diseaseRL is log FC coefficients
			refDB<-diseaseRL
			if(stat=="KS"){refDB<-apply((-1*refDB),2,rank)}
		}
	}else{
		
			refDB=customRefDB
	}
	refDB<-as.matrix(refDB)	
	data<-as.matrix(data)
	#check to make sure the genes in the input profile match those in the reference data set.
	lt<-sum(rownames(data)%in%rownames(refDB))
	if(lt!=nrow(data)) stop('Rownames of input data and reference data do not match')
	#select only those in the input in the ref set too
	refDB<-refDB[rownames(data),]
	if(stat=="KS"){
		#need to have the non rank data first incase want to use the Gant methodology.
		#create two sets of gene names for the input data. One with the highest increase DE expression at the top (rankdata) and one with the largest decrease in DE at the top (rankdatainv)
		datarl<-apply(data,2,function(x) as.integer(rank(x)))
		datarl<-as.matrix(datarl)
		rownames(datarl)<-rownames(data)
		colnames(datarl)<-colnames(data)
	
		rankdata<-.ranknamemat(datarl,2,TRUE)
		rankdatainv<-.ranknamemat(datarl,2,FALSE)
	}else{
		#get the ranks and the sign matrix
		rankdata<-data
		rankdatainv<-data
	}
	
	#type = fixed for a fixed number of genes to use in calculating the KS running sum statistics.
	if(type=="fixed"){
		#calculate enrichment scores, side indicates whether or not to do a one or two tailed test, k is the gene set size.
		outk<-.getESscores(rankdata,rankdatainv,drugPRLs=refDB,uplength=lengthtest,downlength=lengthtest,stat=stat)	
		#get permutation based empirical p-values for enrichment scores	
		signif<-.significantES(ESscores=outk,rankdata=rankdata,rankdatainv=rankdatainv,noperm=noperm,fdr=signif.fdr,uplength=lengthtest,downlength=lengthtest,drugPRLs=refDB,adj=adj,stat=stat)
		names(signif)<-colnames(data)
	}
	#type = range, calculates the KS running sum statistics for a range of gene set sizes.
	if(type=="range"){
		#do dynamic set for find drugs
		outk<-list(length=length(ranges))
		signif<-list(length=length(ranges))
		for(i in 1:length(ranges)){
			#calculate enrichment scores
			outk[[i]]<-.getESscores(rankdata,rankdatainv,drugPRLs=refDB,uplength=ranges[i],downlength=ranges[i],stat=stat)
			#get permutation based empirical p-values for enrichment scores
			signiftemp<-.significantES(ESscores=outk[[i]],rankdata=rankdata,rankdatainv=rankdatainv,noperm=noperm,fdr=signif.fdr,uplength=ranges[i],downlength=ranges[i],drugPRLs=refDB,adj=adj,stat=stat)			
			
			signif[[i]]<-signiftemp	
		}
		
		names(outk)<-paste("range:",ranges,sep="")
		names(signif)<-paste("range:",ranges,sep="")
	}
		
	if(type=="dynamic"){
		#type= dynamic, meaning that the size of the gene set used in the enrichment score calculation from the input profiles, is determined using the adjusted p-values of differential expression
		pvalues<-as.matrix(pvalues)
		#adjust the p-values for multiple hypothesis testing using either Benjamini-Hochberg (BH) or Tibshirani q-value method
    	if(adj=="BH"){
      		pvalues<-apply(pvalues,2,p.adjust)
    	}else{
      		try(pvals<-apply(pvalues,2,function(x) qvalue(x)$qvalues))
      		if(class(pvals)[1]=="try-error"){
      			pvalues<-apply(pvalues,2,p.adjust)
      		}else{
      			pvalues<-pvals	
      		}
    	}
		#signifgenes contains a logical value indicating if have significantly differentially expressed genes after Multiple hypothesis correction for given false discovery rate (dynamic.fdr)
		signifgenes<-apply(pvalues,2,function(x) x<dynamic.fdr)
		#get the ranklists and inverse ranked lists for the reference data set (assuming that the reference set is a rank of decreasing order)
		ranknamedrugPRLs<-.ranknamemat(refDB,2,TRUE)
		ranknamedrugPRLsb<-.ranknamemat(refDB,2,FALSE)
		if(ncol(pvalues)>1){
			#have more than one profile to analyse
			signif<-list(length=ncol(pvalues))
			outk<-list(length=ncol(pvalues))
			positivematrix<-apply(data,2,function(x) x>0)
			upsignif<-positivematrix*signifgenes
			for(i in 1:ncol(pvalues)){
				if(sum(signifgenes[,i])>0){
					#calculate enrichment scores
					
					uplength<-sum(upsignif[,i])
					downlength<-sum(signifgenes[,i])-uplength
					if(uplength<50){
						warning("Number of significantly up-regulated genes is less than 50, top 50 used instead")
						uplength<-50
					}
					if(downlength<50){
						warning("Number of significantly down-regulated genes is less than 50, top 50 used instead")
						downlength<-50
						
					}
					outk[[i]]<-.getESscores(rankdata[,i],rankdatainv[,i],uplength=uplength,downlength=downlength,drugPRLs=refDB,stat=stat)
					#get permutation based empirical p-values for enrichment scores
					signif[[i]]<-.significantES(ESscores=outk[[i]],uplength=uplength,downlength=downlength,rankdata=rankdata,rankdatainv=rankdatainv,noperm=noperm,fdr=signif.fdr,drugPRLs=refDB,adj=adj,stat=stat)}
				else{
					warning(paste("No significantly differentially expressed genes for experiment",colnames(pvalues)[i]))
					outk[[i]]<-NULL
					signif[[i]]<-NULL
				}
			}
			names(outk)<-colnames(data)
			names(signif)<-colnames(data)
		}else{
			#one input profile to analyse
			if(sum(signifgenes)>0){
				positivematrix<-data>0
				upsignif<-positivematrix*signifgenes
				
			   uplength<-sum(upsignif)
			   downlength<-sum(signifgenes)-uplength
			   if(uplength<50){
						warning("Number of significantly up-regulated genes is less than 50, top 50 used instead")
						uplength<-50
					}
				if(downlength<50){
						warning("Number of significantly down-regulated genes is less than 50, top 50 used instead")
						downlength<-50
						
				}
				#calculate enrichment scores
				outk<-.getESscores(rankdata,rankdatainv,drugPRLs=refDB,uplength=uplength,downlength=downlength,stat=stat)	
				#get permutation based empirical p-values for enrichment scores		
				signif<-.significantES(ESscores=outk,rankdata=rankdata,rankdatainv=rankdatainv,noperm=noperm,fdr=signif.fdr,uplength=uplength,downlength=downlength,drugPRLs=refDB,adj=adj,stat=stat)
			}else{
				stop(paste("No significantly differentially expressed genes, consider changing the FDR threshold."))
			}
		names(signif)<-colnames(data)
		names(outk)<-colnames(data)
		}		
	}	
	#so significance (signif) is a list with length=number of input profiles (can be a nested list if type is range so more than one result per profile)

	return(list(scores=outk,significance=signif))	
}
.significantES <-
function(ESscores,rankdata=NULL,rankdatainv=NULL,uplength=250,downlength=250,noperm=1000,fdr=0.01,drugPRLs,adj=adj,stat=c("KS","WSR")){

	#find the significance of the enrichment scores calculated by comparing to random (null hypothesis) permutations of the reference data base

	stat<-match.arg(stat)


	if(is.list(ESscores)){
		#then we have a range of gene set sizes with corresponding enrichment scores	
		resultslist<-list(length=length(ESscores))
		for(i in 1:length(ESscores)){
		#ESscores is a matrix with rows= input profiles and columns = reference profiles
			scorematrix<-ESscores[[i]]
			
				if(is.matrix(scorematrix)){
			
					signifcompounds<-list(length=nrow(scorematrix))
					for(i in 1:nrow(scorematrix)){
      					signifcompounds[[i]]<-.findSignifCompounds(scorematrix[i,],rankdata=rankdata[,i],rankdatainv=rankdatainv[,i],uplength=uplength,downlength=downlength,noperm=noperm,fdr=fdr,adj=adj,stat=stat)
					}
					names(signifcompounds)<-rownames(scorematrix)
				}else{
					signifcompounds<-.findSignifCompounds(scorematrix,rankdata=rankdata,rankdatainv=rankdatainv,uplength=uplength,downlength=downlength,noperm=noperm,fdr=fdr,adj=adj,stat=stat)
				}
			resultslist[[i]]<-signifcompounds
		}
		names(resultslist)<-names(ESscores)
	}else{
		scorematrix<-ESscores
		
		if(is.matrix(scorematrix)){
			#have the diseases in the rows and the drugs in the columns
			#have more than one disease profile:
			signifcompounds<-list(length=nrow(scorematrix))
			for(i in 1:nrow(scorematrix)){
				signifcompounds[[i]]<-.findSignifCompounds(scorematrix[i,],rankdata=rankdata[,i],rankdatainv=rankdatainv[,i],uplength=uplength,downlength=downlength,noperm=noperm,fdr=fdr,adj=adj,stat=stat)
			}		
			names(signifcompounds)<-rownames(scorematrix)
		}else{
			signifcompounds<-list(length=1)
			signifcompounds[[1]]<-.findSignifCompounds(scorematrix,rankdata=rankdata,rankdatainv=rankdatainv,uplength=uplength,downlength=downlength,noperm=noperm,fdr=fdr,adj=adj,stat=stat)
			names(signifcompounds)<-names(ESscores)
		}
	
	resultslist<-signifcompounds
	
}
	#this always outputs a list
	return(resultslist)
	
}
.findSignifCompounds <-
function(scores,rankdata=NULL,rankdatainv=NULL,noperm=100,uplength=250,downlength=250,fdr=fdr,adj=adj,stat=c("KS","WSR")){

	

	if(is.null(rankdata)){stop("Rank data needed for permutation tests")}
		rankdata<-as.matrix(rankdata)
		#generate a random set of gene expression profiles
		size<-nrow(rankdata)
		smatrix<-sapply(1:noperm,function(x) sample.int(size))
		rownames(smatrix)<-rankdata[,1]
		#create the ranklists and decreasing order ranked lists
		samplematrix<-.ranknamemat(smatrix,2,TRUE)
		samplematrixinv<-.ranknamemat(smatrix,2,FALSE)
	

	if(stat=="KS"){
	
		PermScoreup<-.KS.drug.disease(refmatrix=samplematrix,inputmatrix=rankdata,lengthtest=uplength)
		PscoreUpref<-.KS.drug.disease(refmatrix=rankdata,inputmatrix=samplematrix,lengthtest=uplength)
		PscoreUpref<-t(PscoreUpref)
		PermScoredown<-.KS.drug.disease(refmatrix=samplematrix,inputmatrix=rankdatainv,lengthtest=downlength)
		PscoreDownref<-.KS.drug.disease(refmatrix=rankdata,inputmatrix=samplematrixinv,lengthtest=downlength)
		PscoreDownref<-t(PscoreDownref)
	
		#PscoreUp<-(PermScoreup$input_ref+PermScoreup$input_ref)/2
		#PscoreDown<-(PermScoredown$ref_input+PermScoredown$ref_input)/2
		
	
		
		score1<-(PermScoreup-PermScoredown)/2
		score2<-(PscoreUpref-PscoreDownref)/2
		
		check1<-sign(PermScoreup)==sign(PermScoredown)
		check2<-sign(PscoreUpref)==sign(PscoreDownref)
		
		score1[check1]<-0
		score2[check2]<-0		
		
		#PscoreUp<-(PermScoreup$input_ref+PermScoreup$input_ref)/2
		#PscoreDown<-(PermScoredown$ref_input+PermScoredown$ref_input)/2
		
		Pscores<-score1+score2
		Pscores<-Pscores/2
		
		
		#if(combine=="mean"){
				
		#	Pscore<-(PscoreUp+PscoreDown)/2
			
		#}else{
			#combine scores by taking the maximum
			#m<-PscoreUp>PscoreDown
			#Pscore<-PscoreDown
			#Pscore[m]<-PscoreUp[m]
		#	m<-PscoreUp>PscoreDown
		#	Pscore<-PscoreDown
		#	Pscore[m]<-PscoreUp[m]
		#}
	
	}else{
		#using the WSR score (do we want to rank the absolute values?)
		Samplesigns<-sapply(1:noperm,function(x) sample(c(-1,1),size=size,replace=TRUE))
		samplemat<-smatrix*Samplesigns
		rownames(samplemat)<-rownames(rankdata)
		Pscores<-.WSR(refmatrix=samplemat,inputmatrix=rankdata,uplength=uplength,downlength=downlength)
		
		
	}	
		#now check for significance
		#scores are the vector of ES of the input profile(s) versus the reference data set
		#Pscore contains the scores of the input profiles(s) versus the random perturbation data
		#need to separate positive and negative
		#score is always passing one input profile at a time
#ecdfspos<-ecdf(Pscores[which(Pscores>0)])
#		ecdfsneg<-ecdf(Pscores[which(Pscores<0)])
		
		empcdf<-ecdf(abs(Pscores))
		#get the empirical p-value for the actual scores based on the empirical cdf
		emppvals<-vector(length=length(scores))
#     
#        	scorespos<-1-ecdfspos(scores)
#        	scoresneg<-ecdfsneg(scores)
#        emppvals<-scorespos
#        emppvals[which(scores<0)]<-scoresneg[which(scores<0)]	
        
        emppvals<-1-empcdf(abs(scores))
		
		#	
      	if(adj=="BH"){
			pvaladj<-p.adjust(emppvals)
      	}else{
      	
      		qvals<-qvalue(emppvals)
      		if(is.list(qvals)){
      			pvaladj<-qvals$qvalues
      		}else{
      			pvaladj<-p.adjust(emppvals)
      			warning('Q-values error, used BH correction instead')
      		}
			
				
     	}
	
	sel<-which(pvaladj<fdr)
	signifcompounds<-scores[sel]
	names(signifcompounds)<-names(scores)[sel]

	return(signifcompounds)
	
}
.classifysinglelinkage <-
function(resultslist,customClusters=NULL,case=c("disease","drug"),no.signif=10){
	#resultslist output from significantES
	#if length(resultslist)>1 then we have dynamic set sizes for the KS score
	case<-match.arg(case)

	Clusterlist<-list(length=length(resultslist))
		for(i in 1:length(resultslist)){		
			SignifScores<-resultslist[[i]]
			#see if SignifScores is a list - i.e. more than one disease proifle
			if(is.list(SignifScores)){
				Clust<-lapply(SignifScores,function(x) .findCluster(x,customClusters=customClusters,case=case,no.signif=no.signif))
				names(Clust)<-names(SignifScores)
			}else{
				#just have one profile:
				Clust<-.findCluster(SignifScores,customClusters=customClusters,case=case,no.signif=no.signif)
			}
			Clusterlist[[i]]<-Clust		
		
		}
		names(Clusterlist)<-names(resultslist)
		
	
	return(Clusterlist)
}
.classifyaveragelinkage <-
function(ESscores,statistic=c("mean","median"),customClusters=NULL,case=c("disease","drug")){
	statistic<-match.arg(statistic)
	case<-match.arg(case)
	
	if(is.list(ESscores)){
		#then we have done dynamic list of k's
		ClusterAssignments<-list(length=length(ESscores))
		for(i in 1:length(ESscores)){
			#now need to see if have more than one disease profile or not:
			scorematrix<-ESscores[[i]]
			if(is.matrix(scorematrix)){
				avgscores<-list(length=nrow(scorematrix))				
				for(j in 1:nrow(scorematrix)){
				avgscores[[j]]<-.averagecluster(scorematrix[j,],statistic=statistic,customClusters=customClusters,case=case)}
				names(avgscores)<-rownames(scorematrix)
			}else{
				avgscores<-.averagecluster(scorematrix,statistic=statistic,customClusters=customClusters)				
			}
			ClusterAssignments[[i]]<-avgscores
		}
		names(ClusterAssignments)<-names(ESscores)
	}else{
		#fixed k
		scorematrix<-ESscores
		if(is.matrix(scorematrix)){
			avgscores<-list(length=nrow(scorematrix))				
				for(j in 1:nrow(scorematrix)){
					avgscores[[j]]<-.averagecluster(scorematrix[j,],statistic=statistic,customClusters=customClusters,case=case)
				}
			names(avgscores)<-rownames(scorematrix)				
		}else{				
			avgscores<-list(length=1)
			avgscores<-.averagecluster(scorematrix,statistic=statistic,customClusters=customClusters)
			names(avgscores)<-names(scorematrix)
		}
			ClusterAssignments<-avgscores
				
	}
	#This function always returns a list
	return(ClusterAssignments)
}
.findCluster <-
function(SignifScores,customClusters=NULL,case=c("disease","drug"),no.signif=10){
	case<-match.arg(case)
	#check that we have some significant results otherwise return na
	if(length(SignifScores)==0) return(NA)
	if(is.character(SignifScores))return(NA)
	#check if there are more significant results than have been asked to assign to a cluster
	if(length(SignifScores)>no.signif){cat(paste('Number of Significant results greater than',no.signif, 'Using top' ,no.signif, 'hits - consider using average linkage instead'),fill=TRUE)}
	
	#order by most significant resuts (i.e. the highest absolute Enrichment scores) and take the top no.signif to classify
	absscores<-abs(SignifScores)
	SortScores<-sort(absscores,TRUE)
	SScores<-SignifScores[names(SortScores)]
	SignifScores<-SScores[1:no.signif]
	#load default clusters if customClusters is null, according to whether the input profile was a disease or a drug profile.
	if(is.null(customClusters)){
		
		if(case=="disease"){
			data(drugClusters,package="cMap2data")
			Clusters<-drugClusters
		}else{
			data(diseaseClusters,package="DrugVsDiseasedata")
			Clusters<-diseaseClusters
		}
				 
	}else{
		Clusters=customClusters
		cnames<-colnames(Clusters)
		#check the format of the customClusters
		s1<-which(cnames=="Drug")
		s2<-which(cnames=="Disease")
		if(length(s1)+length(s2)==0){stop('customClusters must contain a column with either "Drug" or "Disease" profiles')}
	}
	signifnames<-names(SignifScores)
	#find the cluster for the significant matches, and output to data frame along with the distance between profiles (1-ES), and name of the profile.
	if(case=="disease"){
		sel<-sapply(signifnames,function(x) which(Clusters[,"Drug"]==x))
		sel<-unlist(sel)
		scores<-SignifScores[names(sel)]
		#distscores<-1-abs(scores)
		distscores<-1-abs(scores)
		#are comparing inverted profiles, so positive scores indicate negative correlation and vice versa 
		
		cors<-ifelse(scores>0,-1,1)
		results.clust<-data.frame(names(sel),distscores,Clusters[sel,"Cluster"],cors)
		colnames(results.clust)<-c("Drug","ES Distance","Cluster","RPS")
	}else{
		sel<-sapply(signifnames,function(x) which(Clusters[,"Disease"]==x))
		sel<-unlist(sel)
		scores<-SignifScores[names(sel)]
		#distscores<-1-abs(scores)
		distscores<-1-abs(scores)
		
		cors<-ifelse(scores>0,-1,1)
		results.clust<-data.frame(names(sel),distscores,Clusters[sel,"Cluster"],cors)
		colnames(results.clust)<-c("Disease","ES Distance","Cluster","RPS")
	}
		
return(results.clust)
	
}
.averagecluster <-
function(scores,statistic="mean",customClusters=NULL,case=c("disease","drug"),stat){
	case<-match.arg(case)
	#load default clusters if customClusters is null
	if(is.null(customClusters)){
	
		if(case=="disease"){
			data(drugClusters,package="cMap2data")
			Clusters<-drugClusters
			}else{
			data(diseaseClusters,package="DrugVsDiseasedata")
			Clusters<-diseaseClusters
			}
				 
	}else{
		Clusters=customClusters
	}
	#get the number of clusters available
	clusternumbers<-unique(Clusters[,"Cluster"])
	numberclusters<-length(clusternumbers)
	avgscores<-rep(0,numberclusters)
	medscores<-rep(0,numberclusters)
	avgscoresneg<-rep(0,numberclusters)
	medscoresneg<-rep(0,numberclusters)
	scorepos<-scores>0
	#for each cluster calculate the mean and median scores
	for(k in 1:numberclusters){
		if(case=="disease"){
			selc<-Clusters[which(Clusters[,"Cluster"]==clusternumbers[k]),"Drug"]
		}else{
			selc<-Clusters[which(Clusters[,"Cluster"]==clusternumbers[k]),"Disease"]
			}
			p<-names(scorepos)[scorepos]
			n<-names(scorepos)[!scorepos]
			sp<-selc[selc%in%p]
			sn<-selc[selc%in%n]
		if(length(sp)>0){	
			avgscores[k]<-mean(scores[sp])		
			medscores[k]<-median(scores[sp])
		}
		if(length(sn)>0){
			avgscoresneg[k]<-mean(abs(scores[sn]))
			medscoresneg[k]<-median(abs(scores[sn]))				
		}	
	}
	#now assign to the cluster with the maximum score:
	if(statistic=="mean"){
		cluster<-which.max(avgscores)
		clusterneg<-which.max(avgscoresneg)
	}else{
		cluster<-which.max(medscores)
		clusterneg<-which.max(medscoresneg)
	}
	#now return the scores to elements in the max cluster, sel has the names of the drugs
		if(case=="disease"){
			sel<-Clusters[which(Clusters[,"Cluster"]==clusternumbers[cluster]),"Drug"]
			selneg<-Clusters[which(Clusters[,"Cluster"]==clusternumbers[clusterneg]),"Drug"]	
		}else{
			sel<-Clusters[which(Clusters[,"Cluster"]==clusternumbers[cluster]),"Disease"]
			selneg<-Clusters[which(Clusters[,"Cluster"]==clusternumbers[clusterneg]),"Disease"]
		}
	selscores<-scores[as.character(sel)]
	selscoresneg<-scores[as.character(selneg)]
	distscores<-1-abs(selscores)
	distscoresneg<-1-abs(selscoresneg)
	if(stat=="KS"){
	results<-data.frame(as.character(sel),distscores,rep(clusternumbers[cluster],length(sel)),rep(-1,length(sel)))
	resultsneg<-data.frame(as.character(selneg),distscoresneg,rep(clusternumbers[clusterneg],length(selneg)),rep(1,length(selneg)))}else{
	results<-data.frame(as.character(sel),distscores,rep(clusternumbers[cluster],length(sel)),rep(1,length(sel)))
	resultsneg<-data.frame(as.character(selneg),distscoresneg,rep(clusternumbers[clusterneg],length(selneg)),rep(-1,length(selneg)))
		
		}
	
	#give column names according to whether we have a disease or drug input profile
	if(case=="disease"){
		colnames(results)<-c("Drugs","ES Distance","Cluster","RPS")
		colnames(resultsneg)<-c("Drugs","ES Distance","Cluster","RPS")
	}else{
		colnames(results)<-c("Disease","ES Distance","Cluster","RPS")
		colnames(resultsneg)<-c("Disease","ES Distance","Cluster","RPS")
	}
	return(rbind(results,resultsneg))
}
.getESscores <-
function(diseasematrixup,diseasematrixdown,uplength=250,downlength=250,drugPRLs,stat=c("KS","WSR")){
	stat<-match.arg(stat)

	if(stat=="KS"){
	#get the ranklists and inverse ranked lists for the reference data set (assuming that the reference set is a rank of decreasing order)
	#drugPRLs is in rank decreasing order
	drugPRLs<-abs(drugPRLs)
	ranknamedrugPRLs<-.ranknamemat(drugPRLs,2,TRUE)
	ranknamedrugPRLsb<-.ranknamemat(drugPRLs,2,FALSE)
	#do the "top" of the input profiles to the "bottom" of the reference data set
	upprof<-.KS.drug.disease(refmatrix=ranknamedrugPRLs,inputmatrix=diseasematrixup,lengthtest=uplength)
	upref<-.KS.drug.disease(refmatrix=diseasematrixup,inputmatrix=ranknamedrugPRLs,lengthtest=uplength)
	upref<-t(upref)
	#now compare the "bottom" of the input profiles to the "top" of the reference data set
	downprof<-.KS.drug.disease(refmatrix=ranknamedrugPRLs,inputmatrix=diseasematrixdown,lengthtest=downlength)
	downref<-.KS.drug.disease(refmatrix=diseasematrixup,inputmatrix=ranknamedrugPRLsb,lengthtest=downlength)
	downref<-t(downref)
	
	
	score1<-(upprof-downprof)/2
	score2<-(upref-downref)/2
	
	check1<-sign(upprof)==sign(downprof)
	check2<-sign(upref)==sign(downref)
	score1[check1]<-0
	score2[check2]<-0
	scores<-score1+score2
	#ESinput_ref<-(upprof_ref+downprof_ref)/2
	#TESin_ref<-(upprof_ref+downprof_ref)/2
	#then we have compared the gene sets as defined from the reference profiles to the input profiles as well
	
	#upprof_input<-updiseasedowndrug$ref_input
	#downprof_input<-downdiseaseupdrug$ref_input
	#ESref_input<-(upprof_input+downprof_input)/2
	#TESref_in<-(upprof_input+downprof_input)/2
	
	#combine taking either the mean or the maximum of the scores
	#if(combine=="mean"){
		#all<-(ESinput_ref+ESref_input)/2
	#	all<-(TESin_ref+TESref_in)/2
	#}else{	
		#m<-ESref_input>ESinput_ref
		#all<-ESinput_ref
		#all[m]<-ESref_input[m]
	#	m<-TESref_in>TESin_ref
	#	all<-TESin_ref
	#	all[m]<-TESref_in[m]
	#}

	all<-scores/2
	all<-as.matrix(all)
	

	#all will have the disease profiles as rows and the columns as drugs
	#need to change KS.drug.disease to deal with vectors (i.e one disease profile only)		
	}else{
		#going to do the Gant version instead.
		
		
		all<-.WSR(refmatrix=drugPRLs,inputmatrix=diseasematrixup,uplength=uplength,downlength=downlength)
		all<-as.matrix(all)
		
	}
	
	if(is.matrix(diseasematrixup)){
	rownames(all)<-colnames(diseasematrixup)}
	colnames(all)<-colnames(drugPRLs)
	
	return(all)
	
	
	
}

.WSR<-function(refmatrix,inputmatrix,uplength=250,downlength=250){
	#cant use the refmatrix to inputmatrix with this score as we do not have the DE for the cMap data only the ranks
	#assume that refmatrix and inputmatrix contain the ranked signed profiles (e.g. t or coefs)
	#better to split the data into two for up and down regulated, first rank overall:
	refmatrix<-as.matrix(refmatrix)
	inputmatrix<-as.matrix(inputmatrix)
	N<-nrow(refmatrix)
	m<-uplength+downlength
	inputranksignmat<-apply(inputmatrix,2,function(x) rank(abs(x)))
	refranksignmat<-apply(refmatrix,2,function(x) rank(abs(x)))
	nameinput<-.ranknamemat(inputranksignmat,2,TRUE)
	genesets<-nameinput[1:m,]
	genesets<-as.matrix(genesets)
	#ranks for input, should have max value m (not N) i.e split the data
	
	#rankinput<-seq(lengthtest:1)
	#rankref<-apply(refmatrix,2,rank)
	stat<-matrix(nrow=ncol(inputmatrix),ncol=ncol(refmatrix))
	for(i in 1:ncol(inputmatrix)){
		inputranks<-rank(abs(inputmatrix[genesets[,i],i]))
		
		#refvals<-sign(refmatrix[genesets[,i],])*rankref[genesets[,i],]
		#inputvals<-sign(inputmatrix[genesets[,i],i])*rankinput
		
		inputvals<-inputranks*sign(inputmatrix[genesets[,i],i])
		refvals<-refranksignmat[genesets[,i],]*sign(refmatrix[genesets[,i],])
		statvec<-refvals*inputvals
		stat[i,]<-apply(statvec,2,sum)
		

	}

	maxscore<-sum((N-1:m+1)*(m-1:m+1))

	stat<-stat/maxscore
	#*-1 to invert the scores to be consistent with the KS scores, just for reporting the RPS sign at the end. To clean up later.
	return(stat*-1)
	
}
.KS.drug.disease <-
function(refmatrix,inputmatrix,lengthtest=250){
	
	
	refmatrix<-as.matrix(refmatrix)
	inputmatrix<-as.matrix(inputmatrix)
	N<-lengthtest
	Nh<-nrow(refmatrix)
	nototest<-ncol(inputmatrix)
	#outmat<-matrix(nrow=nototest,ncol=ncol(refmatrix))
	#outmat2<-matrix(nrow=nototest,ncol=ncol(refmatrix))
	outmat<-c()
	#outmat2<-c()
	for(i in 1:nototest){
		
		#for a particular disease signature
		toptail<-inputmatrix[1:lengthtest,i]
		#toptoref contains the references of elements in the refmatrix which are in the gene set of the input profile
		toptoref<-refmatrix%in%toptail


		In<-matrix(0,nrow=nrow(refmatrix),ncol=ncol(refmatrix))
		Out<-matrix(0,nrow=nrow(refmatrix),ncol=ncol(refmatrix))
		#assing 1/N to those elements in the gene set and 1/(Nh-N) to those which are not
		In[toptoref]<-1/N
		Out[!toptoref]<-1/(Nh-N)
		#calculate the cumulative (running sums)
		if(ncol(refmatrix)==1){
				In<-as.matrix(In,ncol=1)
				Out<-as.matrix(Out,ncol=1)
				
		}
		
		In<-apply(In,2,cumsum)
		Out<-apply(Out,2,cumsum)
		#take the difference to get the running Enrichment scores
		Testmat<-In-Out
		#take the maxium as the final enrichment score
		refvalues<-apply(Testmat,2,function(x) x[which.max(abs(x))])
		#outmat[i,]<-apply(Testmat,2,max)
		#outmat[i,]<-refvalues
		outmat<-rbind(outmat,refvalues)
	}
	
		#if we are also comparing the gene sets from the reference data to the input profiles
		#for(j in 1:ncol(refmatrix)){
		
		#	toptail<-refmatrix[1:lengthtest,j]
		#	toptoinput<-inputmatrix%in%toptail
		#	In<-matrix(0,nrow=nrow(refmatrix),ncol=ncol(inputmatrix))
		#	Out<-matrix(0,nrow=nrow(refmatrix),ncol=ncol(inputmatrix))

		#	In[toptoinput]<-1/N
		#	Out[!toptoinput]<-1/(Nh-N)
			
			#in case there is only one reference profile being classified
		#	if(ncol(inputmatrix)==1){
		#		In<-as.matrix(In,ncol=1)
		#		Out<-as.matrix(Out,ncol=1)
				
		#	}
		#	In<-apply(In,2,cumsum)
		#	Out<-apply(Out,2,cumsum)
		#	Testmat<-In-Out
			#outmat2[,j]<-apply(Testmat,2,max)
		#	refvalues2<-apply(Testmat,2,function(x) x[which.max(abs(x))])
			#outmat2[,j]<-refvalues2
		#	outmat2<-cbind(outmat2,refvalues2)
		#}
		#rownames(outmat2)<-colnames(inputmatrix)
		#colnames(outmat2)<-colnames(refmatrix)

	
	
	rownames(outmat)<-colnames(inputmatrix)
	colnames(outmat)<-colnames(refmatrix)
	
	#return(list(input_ref=outmat,ref_input=outmat2))
	return(outmat)
}
.ranknamemat <-
function(matrix,ref,decreasing){
		#rank the matrix, returning the names, starting with the most signif
		
			#want to start with largest first (e.g. rank expression values) 
			if(is.matrix(matrix)){
			rankedmat<-apply(matrix,ref,function(x) names(x)[sort.list(as.integer(x),decreasing=decreasing,method="radix")])
			}else{
				rankedmat<-names(sort(matrix,decreasing=decreasing))
			}
		#output is a vector or matrix of ranked lists sorted according to parameter 'decreasing' and with rank values replaced by their names
		return(rankedmat)
		
}
.writecytoscape <-
function(case=c("disease","drug"),clusterassignments,filename,customsif=NULL,customedge=NULL,nodename=NULL,customClusters=NULL){
	
	case<-match.arg(case)


	#load default networks and generate SIF and Edge attribute files if custom networks have not been provided
	if(!is.null(customsif)){
		networksif=customsif
		cytoedge=customedge
		RPS<-cytoedge
		RPS[,5]<-1
		if(case=="disease"){
			nodetype<-data.frame(cbind(unique(as.character(RPS[,1])),rep("=",length(unique(RPS[,1]))),rep("drug",length(unique(RPS[,1])))))
			
		}else{
			nodetype<-data.frame(cbind(unique(as.character(RPS[,1])),rep("=",length(unique(RPS[,1]))),rep("disease",length(unique(RPS[,1])))))
			}
			
		Clusters<-customClusters
		labels<-NULL
	}else{
		if(case=="disease"){
			data(cytodrug,package="cMap2data")
			data(drugClusters,package="cMap2data")
			Clusters<-drugClusters
			cytooutput<-cytodrug
			nodetype<-data.frame(cbind(as.character(Clusters[,1]),rep("=",nrow(Clusters)),rep("drug",nrow(Clusters))))
			data(druglabels,package="cMap2data")
			labels<-druglabels
		}else{
			data(cytodisease,package="DrugVsDiseasedata")
			data(diseaseClusters,package="DrugVsDiseasedata")
			Clusters<-diseaseClusters
			nodetype<-data.frame(cbind(as.character(Clusters[,1]),rep("=",nrow(Clusters)),rep("disease",nrow(Clusters))))
			cytooutput<-cytodisease
			data(diseaselabels,package="DrugVsDiseasedata")
			labels<-diseaselabels
			
		}
		#create SIF format of nodename1,1,nodename2. N.B. the edge type is always 1
		suppressWarnings(networksif<-data.frame(cbind(as.character(cytooutput[,1]),rep(1,nrow(cytooutput)),as.character(cytooutput[,2]))))
		#create edge attribute file: nodename1, "(1)", nodename2, "=", Distance between node1 and node2
		suppressWarnings(cytoedge<-data.frame(cbind(as.character(cytooutput[,1]),rep("(1)",nrow(cytooutput)),as.character(cytooutput[,2]),rep("=",nrow(cytooutput)),as.character(cytooutput[,3]))))
		suppressWarnings(RPS<-data.frame(cbind(as.character(cytooutput[,1]),rep("(1)",nrow(cytooutput)),as.character(cytooutput[,2]),rep("=",nrow(cytooutput)),as.character(cytooutput[,4]))))
	
		
	}
	cytonode<-as.data.frame(cbind(as.character(Clusters[,1]),rep("=",nrow(Clusters)),Clusters[,2]))	
	#create the sif, edge files for the reference network
	edgeheader<-vector(mode="character",length=5)
	edgeheader[1]<-"Distance"
	nodeheader<-vector(mode="character",length=3)
	nodeheader[1]<-"ClusterID"
	#RPS = Running-sum Peak Sign
	nodetypeheader<-vector(mode="character",length=3)
	nodetypeheader[1]<-"NodeType"
	RPSheader<-vector(mode="character",length=5)
	RPSheader[1]<-"RPS"
	write.table(networksif,file=paste(filename,"network.sif",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=FALSE)
	write.table(cytoedge,file=paste(filename,"Distance.EA",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=edgeheader)
	write.table(cytonode,file=paste(filename,"Cluster.noa",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=nodeheader)
	write.table(RPS,file=paste(filename,"RPS.EA",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=RPSheader)
	write.table(nodetype,file=paste(filename,"NodeType.noa",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=nodetypeheader)
	
	signifmatches<-clusterassignments
	#for the significant matches and genereate the sif and edge information
	allnames<-c()
	inputsif<-c()
	inputedge<-c()
	inputRPS<-c()
	for(i in 1:length(signifmatches)){
		res<-signifmatches[[i]]
		if(!is.data.frame(res)){

			for(j in seq(length=length(res))){
				t<-as.matrix(res[[j]])
				if(!is.null(t)){
					nodename=paste(names(signifmatches)[i],names(res)[[j]],sep="_")
					m<-data.frame(cbind(rep(nodename,nrow(t)),as.character(t[,1])))
					inputsiftemp<-data.frame(cbind(as.character(m[,1]),rep(1,nrow(m)),as.character(m[,2])))	
					inputedgetemp<-data.frame(cbind(as.character(m[,1]),rep("(1)",nrow(m)),as.character(m[,2]),rep("=",nrow(m)),t[,2]))
					inputRPStemp<-data.frame(cbind(as.character(m[,1]),rep("(1)",nrow(m)),as.character(m[,2]),rep("=",nrow(m)),t[,4]))
					colnames(inputsiftemp)<-colnames(networksif)
					colnames(inputedgetemp)<-colnames(cytoedge)
					colnames(inputRPStemp)<-colnames(cytoedge)
					allnames<-c(allnames,nodename)
				}
				inputsif<-rbind(inputsif,inputsiftemp)
				inputedge<-rbind(inputedge,inputedgetemp)
				inputRPS<-rbind(inputRPS,inputRPStemp)
			}
			
		}else{
			nodenames<-names(signifmatches)
			allnames<-c(allnames,nodenames)
			if(!is.null(res)){
				res<-as.matrix(res)
				m<-data.frame(cbind(rep(nodename[i],nrow(res)),as.character(res[,1])))
				inputsiftemp<-data.frame(cbind(as.character(m[,1]),rep(1,nrow(m)),as.character(m[,2])))
				inputedgetemp<-data.frame(cbind(as.character(m[,1]),rep("(1)",nrow(m)),as.character(m[,2]),rep("=",nrow(m)),res[,2]))
				inputRPStemp<-data.frame(cbind(as.character(m[,1]),rep("(1)",nrow(m)),as.character(m[,2]),rep("=",nrow(m)),res[,4]))
				colnames(inputsiftemp)<-colnames(networksif)
				colnames(inputedgetemp)<-colnames(cytoedge)
				colnames(inputRPStemp)<-colnames(cytoedge)

			}
		}
		inputsif<-rbind(inputsif,inputsiftemp)
		inputedge<-rbind(inputedge,inputedgetemp)
		inputRPS<-rbind(inputRPS,inputRPStemp)
		
	}
	#cluster for new nodes:
	inputnode<-data.frame(cbind(unique(as.character(allnames)),rep("=",length(unique(allnames))),rep(0,length(unique(allnames)))))
	inputNT<-data.frame(cbind(unique(as.character(allnames)),rep("=",length(unique(allnames))),rep(case,length(unique(allnames)))))
	colnames(inputnode)<-colnames(cytonode)
	#append edges and attributes of the input profiles to nodes in the reference set file.
	write.table(inputsif,file=paste(filename,"network.sif",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=FALSE,append=TRUE)
	write.table(inputedge,file=paste(filename,"Distance.EA",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=FALSE,append=TRUE)
	write.table(inputnode,file=paste(filename,"Cluster.noa",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=FALSE,append=TRUE)
	write.table(inputRPS,file=paste(filename,"RPS.EA",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=FALSE,append=TRUE)
	write.table(inputNT,file=paste(filename,"NodeType.noa",sep="_"),quote=FALSE,na="",row.names=FALSE,col.names=FALSE,append=TRUE)
	if(!is.null(labels)){
		#now do the mapping for the drug/disease urls
		labels<-data.frame(cbind(as.character(labels[which(labels[,1]%in%Clusters[,1]),1]),rep("=",length(unique(allnames))),as.character(labels[which(labels[,1]%in%Clusters[,1]),2])))
		nulllabels<-data.frame(cbind(unique(as.character(allnames)),rep("=",length(unique(allnames))),rep(NA,length(unique(allnames)))))
		colnames(nulllabels)<-colnames(labels)
		labelout<-rbind(labels,nulllabels)
		lname<-vector(mode="character",length=3)
		lname[1]<-"NodeURL"
		write.table(labelout,file=paste(filename,"NodeURL.noa",sep="_"),quote=FALSE,row.names=FALSE,col.names=lname)
		
	}
	
}
.combineEnsembl <-
function(data,annotation,genelist,type=c("average","medpolish","maxvar")){
  	#remove probes which map to more than one gene...
	countprobe<-table(annotation[,1])
	selp<-names(which(countprobe==1))
	data<-as.matrix(data[selp,])
	annotation<-annotation[which(annotation[,1]%in%selp),]
	#find out how many probes map to a given gene
	countgene<-table(annotation[,2])
	#find all those which are 1-1 probe-gene mappings
	sel<-names(which(countgene==1))
	probes<-annotation[which(annotation[,2]%in%sel),1]
	temp<-as.matrix(data[probes,])
	rownames(temp)<-annotation[which(annotation[,2]%in%sel),2]
	#find those which are not a 1-1 mapping
	combinegenes<-names(which(countgene!=1))
	subdata<-as.matrix(data[-(which(rownames(data)%in%probes)),])
	subannot<-annotation[which(annotation[,1]%in%rownames(subdata)),]
	#sort data according to the genes they map to
	sorder<-subannot[order(subannot[,2]),]
	sdorder<-subdata[sorder[,1],]
  	sdorder<-as.matrix(sdorder)
	tabann<-table(sorder[,2])
	#combine results where there is not a 1-1 mapping
	if(type=="average"){
		other<-.combined(tabann,sdorder)}
	if(type=="medpolish"){
		other<-.combinemed(tabann,sdorder)
	}
	if(type=="maxvar"){
		cumulativecount<-cumsum(tabann)
		other<-.maxvar(sdorder,cumulativecount)
	}
	if(type=="max"){
		other<-convertProbes(subdata,subannot)
	}
	rownames(other)<-combinegenes
	d<-rbind(temp,other)
	l1<-length(which(rownames(d)%in%genelist))
	if(l1!=length(genelist)){cat('Note: Ensembl genes do not match genelist in reference data. Consider uploading pre-processed lists to classifyprofiles',fill=TRUE)}
	d[which(rownames(d)%in%genelist),]
}
.combined <-
function(tabann,sdorder){
	
	refs<-cumsum(tabann)
	other2<-matrix(nrow=(length(tabann)-1),ncol=ncol(sdorder))
  	if(ncol(sdorder)==1){
  	  other<-mean(sdorder[1:refs[1],])
      for(i in 2:length(tabann)){
    	  other2[i-1,]<-mean(sdorder[refs[i-1]:refs[i],]) 
      }
   	}else{
   	  other<-colMeans(sdorder[1:refs[1],])
  	  for(i in 2:length(tabann)){
		  other2[i-1,]<-colMeans(sdorder[refs[i-1]:refs[i],])
		
	  }
   }
	colnames(other2)<-colnames(sdorder)
	
	rbind(other,other2)
}
.combinemed<-
function(tabann,sdorder){
	
	refs<-cumsum(tabann)
	other2<-matrix(nrow=(length(tabann)-1),ncol=ncol(sdorder))
  	if(ncol(sdorder)==1){
  	  mfit<-medpolish(sdorder[1:refs[1],])
  	  other<-mfit$overall+mfit$col
      for(i in 2:length(tabann)){
      	  mfit<-medpolish(sdorder[refs[i-1]:refs[i],])
    	  other2[i-1,]<-mfit$overall+mfit$col 
      }
   	}else{
   		mfit<-medpolish(sdorder[1:refs[1],])
   	  other<-mfit$overall+mfit$col
  	  for(i in 2:length(tabann)){
  	  	  mfit<-medpolish(sdorder[refs[i-1]:refs[i],])
		  other2[i-1,]<-mfit$overall +mfit$col
		
	  }
   }
	colnames(other2)<-colnames(sdorder)
	
	rbind(other,other2)
}

.matrixVar <- function(x) {
  x <- as.matrix(x)
  n <- ncol(x)
  v <- (((x * x) %*% matrix(rep(1, n), nrow=n)) - n * rowMeans(x) * rowMeans(x))/(n - 1)
  v[v < 0] <- 0
  v
}
.maxvar<-function(probevals,cumulativecount){
	
	maxvar<-vector(length=length(cumulativecount))
	maxvar[1]<-rownames(probevals)[which.max(.matrixVar(probevals[1:cumulativecount[1],]))]
	for(i in 2:length(cumulativecount)){
		rnames<-rownames(probevals)[(cumulativecount[i-1]+1):cumulativecount[i]]
		maxvar[i]<-rnames[which.max(.matrixVar(probevals[(cumulativecount[i-1]+1):cumulativecount[i],]))]
	}
	return(probevals[maxvar,])
}
convertProbes<-function(subdata,subannot){
	

	#sort data according to the genes they map to
	#sorder<-subannot[order(subannot[,2]),]
	#sdorder<-subdata[sorder[,1],]
  	#sdorder<-as.matrix(sdorder)
  	annot<-subannot[,2]
  	nexps<-ncol(subdata)
  	noprobes<-nrow(subdata)
  	check1<-subdata[,1:nexps]<(noprobes/2)
	data1<-subdata[,1:nexps]-(noprobes/2)
	data2<-(noprobes/2)-subdata[,1:nexps]
	outranks<-data1[,1:nexps]
	outranks[check1]<-data2[check1]
  	
	out<-apply(outranks,2,function(x) .ptog(x,annot))
		
}

.ptog<-function(data,annot){
	d<-data.frame(cbind(data,annot))
	colnames(d)<-c("data","gene")
	v<-d$gene[i<-order(d$data)]
	ind<-!duplicated(v,fromLast=TRUE)
	idx<-setNames(seq_len(nrow(d))[i][ind],v[ind])
	data<-data[idx]
	names(data)<-names(idx)
	data<-data[order(names(idx))]
	
}




.convertEnsembl <-
function(annotation,values){
	#values should be the list of probes for the selected platform
	sel<-which(annotationlist[,1]==annotation)
	if(length(sel)==0){stop('Annotation file not available for this platform. Must be one of Affymetrix, HGU133A, HGU133A2 or HGU133Plus2. Alternatively upload an annotation file to process the data')}
	annref<-annotationlist[sel,2]
	human_mart<-useMart("ensembl","hsapiens_gene_ensembl")
	annotation_ensembl<-getBM(attributes=c(annref,"hgnc_symbol"),filters=annref,values=values,mart=human_mart)
	#this will be a matrix each attributes are columns rows are values
	return(cbind(annotation_ensembl[[1]],annotation_ensembl[[2]]))
	
}
.fitlms <-
function(normalised,experimentAEpdata=NULL,customfactors=NULL,experiment=NULL,case=c("disease", "drug"),statistic=c("coef", "t","diff")){
		
		case<-match.arg(case)
		statistic<-match.arg(statistic)
		
		#if statistic chosen is diff, simply remove logarithms and use the standard linear model coefficients
		if(statistic=="diff"){
			normalised<-2^(normalised)
			statistic<-"coef"
		}
		
					
			if(is.null(customfactors)){
				if(is.null(experimentAEpdata)){stop("No factor data provided in either experimentAEpdata or customfactors")}
				#find available factors from the meta data
				factorvals<-grep("Factor.Value",ignore.case=TRUE,colnames(experimentAEpdata),value=T)
				if(length(factorvals)==0){factorvals<-grep("FactorValue",ignore.case=TRUE,colnames(experimentAEpdata),value=T)}
				if(length(factorvals)==0){stop("No Factor Values found in the pData. Check the Array Express entry and consider using the customfactors to supply factor values.")}
				if(case=="disease"){
					factordisease<-grep("state",ignore.case=TRUE,factorvals,value=T)
					if(length(factordisease)==0)factordisease<-grep("status",ignore.case=TRUE,factorvals,value=T)
				}else{
						factordisease<-grep("compound",ignore.case=TRUE,factorvals,value=T)
				}
				
				
				if(length(factordisease)==0){
					fv<-c()
					for(i in 1:length(factorvals)){
						fv<-paste(fv,factorvals[i])	
					}
					stop(paste("No disease factor present for experiment, use customfactors to provide one. Available factor values are:",fv))
				}
					factors<-experimentAEpdata[,factorvals]
					numberfactors<-length(factorvals)
					factor<-experimentAEpdata[,factordisease]
					factor<-as.vector(factor)
			
			}else{
				#using customfactors from local made file		
				#if option passed is a character string, assume this is a file path:
				if(is.character(customfactors))customfactors<-read.table(customfactors,header=TRUE)
				#if it is an R data object then extract factor values based on column names:
				if(is.matrix(customfactors)){
					factorvals<-colnames(customfactors)
				}else{
					factorvals<-names(customfactors)
				}
				factors<-customfactors
				#find the "main" explanatory factor
				if(case=="disease"){
					factordisease<-grep("disease",ignore.case=TRUE,factorvals,value=T)
				}else{
					factordisease<-grep("compound",ignore.case=TRUE,factorvals,value=T)
				}
				
				factorfiles<-grep("CEL",ignore.case=TRUE,factorvals,value=T)
				if(length(factorfiles)==0){stop("No CEL file factor present for imported data, check custom factors")}
				if(length(factordisease)==0){stop("No disease/compound factor present for experiment, check custom factors")}
				if(length(factordisease)>1){stop("More than one disease/compound factor found, check custom factors")}
				if(length(factorfiles)>1){warning("More than one CEL file name factor found")}
				#make sure the factors match the order the files will have been read in from the directory by ReadAffy:
				
				factors<-customfactors[order(customfactors[,factorfiles]),]
				factor<-as.vector(factors[,factordisease])
				factors<-factors[,-which(factorvals==factorfiles)]
				factorvals<-factorvals[-which(factorvals=="CEL")]
				numberfactors<-length(factorvals)
				
			}
			
			if(length(unique(factor))<2){
				fv<-c()
				for(i in 1:length(factorvals)){
					fv<-paste(fv,factorvals[i])
				}
	
stop(paste("Disease/Compound factor has less than 2 levels, available factors are:",fv))}
			
			if(!is.factor(factor)){factor<-factor(factor)}
			#generate the main explanatory factore design matrix
			suppressWarnings(designdisease<-model.matrix(~0+factor))
			colnames(designdisease)<-make.names(colnames(designdisease))
			nlevdisease<-nlevels(as.factor(factor))

			design<-c()
			cnames<-c()
			ranklist<-c()
			pvalues<-c()
			factstoinc<-0
		#stratify the main factor according to additional factors present and check if there are enough replicates at each level to fit a linear model
		if(numberfactors>1){
			fref<-which(factorvals==factordisease)
			factnodisease<-as.matrix(factors[,-fref])
			colnames(factnodisease)<-colnames(factors)[-fref]
			factlevstoinc<-vector("list",length=(numberfactors-1))
			factornames<-colnames(factnodisease)
			for(i in 1:(numberfactors-1)){
				
				if(numberfactors>2){
					factorcount<-table(factor,factnodisease[,i])
				}else{
					factorcount<-table(factor,factnodisease)
				}
				minreps<-apply(factorcount,2,min)
				#check have at least 2 replicates
				tu<-which(minreps>1)
				if(length(tu)>0){
					
					factstoinc<-factstoinc+1
					factlevstoinc[[factstoinc]]<-colnames(factorcount)[tu]
					if(is.matrix(factnodisease))factornames[factstoinc]<-colnames(factnodisease)[i]	
				}
				
			}
			names(factlevstoinc)<-factornames
		}
		

		#if there are additional factors to stratify the data by:
			if(factstoinc>0){
			
				
					for(i in seq(length=factstoinc,by=1)){
						if(i>1){
							#store the names of the contrasts
							ncur<-colnames(ranklist)
						}
						#extract the additional factor to be used in the linear model
						uniquefac<-factlevstoinc[[i]]	
	
						nlev<-length(uniquefac)
						cnames<-c()
						for(j in 1:nlev){
							#for each level of the second factor
							if(is.double(normalised)){vec<-vector(length=ncol(normalised))}else{
							vec<-vector(length=dim(normalised)["Samples"])}
							if(is.matrix(factnodisease)){
							#create a vector to store whether or not a replicate (sample) is to be included
							vec[which(factnodisease[,names(factlevstoinc)[i]]==uniquefac[j])]<-1}else{vec[which(factnodisease==uniquefac[j])]<-1}
							#create the sub design matrix
							designtouse<-designdisease*vec
							#fit the linear model
							lmodel<-lmFit(normalised,designtouse)
						
							#fit all available contrasts according to the number of levels of the main factor
							if(nlevdisease>2){
								
								rl<-try(.multcontrast(nlevdisease,designtouse,lmodel,expname=uniquefac[j],statistic=statistic))
								if(class(rl)[1]=="try-error"){warning(paste("Linear model failed for factor",uniquefac[j]))}
								#cnames<-c(cnames,colnames(rl$ranks))
								ranklist<-cbind(ranklist,rl$ranks)
								pvalues<-cbind(pvalues,rl$pvalues)
							}
							else{
								#try making contrast matrix by hand instead and search for "normal" keyword in disease factor
								
								rl<-try(.singlecontrast(normalised,designtouse,expname=uniquefac[j],statistic=statistic))
								if(class(rl)[1]=="try-error"){warning(paste("Linear model failed for factor",uniquefac[j]))}
								#need to do the other outputs e.g. SAM, Q-values, p-values
								ranklist<-cbind(ranklist,rl$ranks)
								pvalues<-cbind(pvalues,rl$pvalues)
								
							}
										
						}	
						
					}
					#now do for main factor only
					rl<-.treatmentonlyfit(nlevdisease,normalised,designdisease,factor,expname="all",statistic=statistic)
					
					ncur<-colnames(ranklist)
					ranklist<-cbind(ranklist,rl$ranks)
					pvalues<-cbind(pvalues,rl$pvalues)
					
					colnames(ranklist)<-c(ncur,colnames(rl$ranks))
					colnames(pvalues)<-c(ncur,colnames(rl$ranks))
								
				}else{
				#have no subtypes to do only the main factor				
					rl<-.treatmentonlyfit(nlevdisease,normalised,designdisease,factor,expname="all",statistic=statistic)		
					ranklist<-rl$ranks
					pvalues<-rl$pvalues			
				}
		
		
		return(list(ranklist=ranklist,pvalues=pvalues))
}
.treatmentonlyfit <-
function(nlevdisease,normalised,designdisease,factor,expname,statistic){
	
	if(nlevdisease>2){
		#fit multiple level contrasts.
		lmodel<-lmFit(normalised,designdisease)
		rl<-try(.multcontrast(nlevdisease,factor,lmodel,expname=expname,statistic=statistic))
		if(class(rl)[1]=="try-error")stop("Unable to fit the linear model, check data inputs")
		rt<-rl$ranks
		pvalues<-rl$pvalues

								
		}else{
			#assume have case/control
			ranklist<-try(.singlecontrast(normalised,designdisease,expname=expname,statistic=statistic))
			if(class(ranklist)[1]=="try-error")stop("Unable to fit the linear model, check data inputs")
			rt<-ranklist$ranks
			pvalues<-ranklist$pvalues
		
		}

	return(list(ranks=rt,pvalues=pvalues))	
	
}
.singlecontrast <-
function(normalised,designdisease,expname,statistic){
		#fit model of DE between case and control for the main factor
		fit<-lmFit(normalised,designdisease)
		#do both possible contrasts (case-control and control-case)
		contrdisease<-c(1,-1)
		contrdisease<-rbind(contrdisease,c(-1,1))
		rownames(contrdisease)<-colnames(designdisease)
		fit2<-contrasts.fit(fit,contrdisease)
		coefs<-fit2$coefficients
		#do empirical bayes procedure from limma, and extract either logFC (coefficients) or t-statistic (t)
		fite<-eBayes(fit2)
		if(statistic=="coef"){
			ranklist<-fite$coefficients
		}else{
			ranklist<-fite$t
		}
		#get the corresponding pvalues
		pvalues<-fite$p.value
		#set the names for each of the contrasts
		cnames<-c(paste(colnames(designdisease)[1],"-",colnames(designdisease)[2],":",expname,sep=""),paste(colnames(designdisease)[2],"-",colnames(designdisease)[1],":",expname,sep=""))
		colnames(ranklist)<-cnames
		colnames(pvalues)<-cnames
				
		return(list(ranks=ranklist,pvalues=pvalues))
	
}
.multcontrast <-
function(nlevdisease,factor,lmodel,expname,statistic){
		#find all combinations of two levels of the main factor
		nocombs<-choose(nlevdisease,2)
		contrmat<-matrix(0,nrow=nocombs,ncol=nlevdisease)
		colnames(contrmat)<-colnames(lmodel$coefficients)
		combinations<-combn(colnames(lmodel$coefficients),2)
		#generate the contrasts for each possible combination
		for(k in 1:nocombs){
			contrmat[k,combinations[1,k]]<-1
			contrmat[k,combinations[2,k]]<-(-1)
		}
		
		fullcontrmat<-rbind(contrmat,contrmat*(-1))
		fullcontrmat<-t(fullcontrmat)
		#fit the contrasts and do empirical bayes from limma
		fitfull<-contrasts.fit(lmodel,fullcontrmat)
		modele<-eBayes(fitfull)
		#extract either log FC (coefficient) or the t-statistic (t)
		if(statistic=="coef"){
			ranklist<-modele$coefficients
		}else{
			ranklist<-modele$t
		}
		#get the pvalues
		pvalues<-modele$p.value
		cnames1<-paste(combinations[1,],"-",combinations[2,],":",expname,sep="")
		cnames2<-paste(combinations[2,],"-",combinations[1,],":",expname,sep="")
		colnames(ranklist)<-c(cnames1,cnames2)
		colnames(pvalues)<-c(cnames1,cnames2)
		return(list(ranks=ranklist,pvalues=pvalues))
									
}

.mrs = function(platform1.data, platform2.data, p1.names=0, p2.names=0, skip.match = FALSE){
######################################################################
#Copyright Jason Rudy & Faramarz Valafar 2009-2010
	
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
	
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
	
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
######################################################################
	
######################################################################
#This file contains code adapted from the WebArrayDB program,
#available from http://www.webarraydb.org/webarray/index.html  Please 
#cite appropriately when using these functions.  Type 
#citation("CONOR") at the R prompt for details.
	
	
#Implements the median rank scores algorithm for cross
#platform normalization.
	

	#This function is basically a wrapper for normalizeGQ
	
	#Match names
	input = .processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
	
	#Prepare for normalizeGQ
	combined = cbind(input$x,input$y)
	pf = c(seq(1,1,length.out=dim(input$x)[2]),seq(2,2,length.out=dim(input$y)[2]))
	
	#Call normalizeGQ
	nmrs = .norm_mrs(combined,pf)
	
	#Split the results and return
	out=split(seq(pf),pf)
	out[[1]] = nmrs[,out[[1]]]
	out[[2]] = nmrs[,out[[2]]]
	names(out) <- c("x","y")
	return(out)
}



.norm_mrs <- function(M, pf, ...) { # median rank  scores algorithm  - a modification of quantile for multi-platform data
	idx <- tapply(seq(pf), pf, function(x) x)
	if (length(pf)<=0) return(M)
	#if (length(idx)==1) return(normalizeQuantiles(M))
	imax <- which.max(sapply(idx, length)) # for reference
	ref_med_srt <- sort(apply(M[, idx[[imax]]], 1, function(x) median(x, na.rm=TRUE) ))
	idx[imax] <- NULL
	lapply(unlist(idx), function(i) M[,i] <<- ref_med_srt[rank(M[,i])])
	invisible(M)
	}
	
	
.processplatforms = function(datalist, namesvec=NULL, skip.match=FALSE){
	#Convert data from various formats to the proper format for use 
	#with all the crossnorm normalization functions
	
	for(i in 1:length(datalist)){
		if(is.matrix(datalist[[i]])){
			datalist[[i]] <- as.data.frame(datalist[[i]])
		}
	}
	
	if (is.null(namesvec)){
		namesvec <- numeric(length(datalist))
		for (i in 1:length(datalist)){
			namesvec[i] <- 0
		}
	}
	
	#Put the row names in their places
	for (i in 1:length(namesvec)){
		if(namesvec[i] != 0){
			rownames(datalist[[i]]) = datalist[[i]][,namesvec[i]]
			datalist[[i]] = datalist[[i]][,-1*namesvec[i],drop=FALSE]
		}	
	}

	if(!skip.match){
		#Create the common genes list
		commongenes <- rownames(datalist[[1]])
		for (i in 2:length(datalist)){
			commongenes <- intersect(commongenes,rownames(datalist[[i]]))
		}
	
	
		#Put it all together
		for (i in 1:length(datalist)){
			datalist[[i]] <- datalist[[i]][commongenes,,drop=FALSE]
		}
	}
	return(datalist)
}

