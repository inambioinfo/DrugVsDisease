generateprofiles <-
function(input=c("AE","GEO","localAE","local"),normalisation=c("rma","mas5"),accession=NULL, customfile=NULL,celfilepath=NULL,sdrfpath=NULL,case=c("disease","drug"),statistic=c("coef", "t","diff"),annotation=NULL,factorvalue=NULL){
	#NB coef as default statistic should be more appropriate as it is more consistent with the generation of the cmap rank profiles which take averages (case, control) and then do ratios.
  #load required packages:

	input<-match.arg(input)
  normalisation<-match.arg(normalisation)
	experimentname=NULL
	statistic<-match.arg(statistic)
	case<-match.arg(case)
	data(genelist,package="DrugVsDiseasedata")
	data(annotationlist,package="DrugVsDiseasedata")
	
	if(input=="GEO"){
		if(is.null(accession)){stop('An accession reference is required to download from GEO database')}
		#download GDS file using GEOquery package
		GEOdownload<-.datafromGEO(accession,factorvalue=factorvalue,case=case)
		#extract the microarray platform
		platform<-GEOdownload$platform
		#check if the platform annotation has been passed as a parameter
		if(!is.null(annotation)){
			platform<-annotation
		}
		if(is.null(platform)){stop('Type of platform unavailable in the experiment information. Please provide using the annotation option')}
		#extact the probenames
		probenames<-rownames(GEOdownload$normalised)
		#find the annotation from BioMart of the probes to Ensembl genes
		annotation<-.convertEnsembl(platform,values=probenames)
		#combine the expression values from probes to genes
		datagenes<-.combineEnsembl(GEOdownload$normalised,annotation,genelist)
		#fit the linear models
		d2<-.fitlms(normalised=datagenes,customfactors=GEOdownload$customfactors,case=case,statistic=statistic,experiment=experimentname)
	    
	}
	if(input=="localAE"){
		if(is.null(celfilepath)|is.null(sdrfpath)){stop('The file path for CEL files and the sdrf file path are required to import Array Express files from a local directory')}
		#read Array Express downloaded CEL and sdrf files, using ArrayExpress package
		AEdownload<-.readlocalAE(celpath=celfilepath,sdrfpath=sdrfpath)
		#normalise data
		normaliseddata<-.normalisedata(AEdownload$affyobject,normalisation=normalisation)
		#extract platform type
		platform<-annotation(normaliseddata)
		#get the expression data
		exprdata<-exprs(normaliseddata)
		if(normalisation=="mas5")exprdata<-log(exprdata,2)
		annotation<-.convertEnsembl(platform,rownames(exprdata))
		datagenes<-.combineEnsembl(exprdata,annotation,genelist)
		d2<-.fitlms(normalised=datagenes,experimentAEpdata=AEdownload$factors,experiment=experimentname,case=case,statistic=statistic)
	}
	if(input=="AE"){
		if(is.null(accession)){stop('An accession reference is required to download from Array Express database')}
		AEdownload<-.datafromAE(experiment=accession)
		platform<-AEdownload$platform

		normaliseddata<-.normalisedata(AEdownload$affyobject,normalisation=normalisation)
		exprdata<-exprs(normaliseddata)
		
		if(normalisation=="mas5")exprdata<-log(exprdata,2)
		annotation<-.convertEnsembl(platform,values=rownames(exprdata))
		datagenes<-.combineEnsembl(exprdata,annotation,genelist)
		d2<-.fitlms(normalised=datagenes,experimentAEpdata=AEdownload$factors,experiment=experimentname,case=case,statistic=statistic)
			
	}
	if(input=="local"){
		if(is.null(celfilepath)|is.null(customfile)){stop('Please specify a path to read the CEL files from using celfilepath and a factor file using customfile')}
		affyobj<-.readlocalCEL(celfilepath)
		normaliseddata<-.normalisedata(affyobj,normalisation=normalisation)
		platform<-annotation(normaliseddata)
		exprdata<-exprs(normaliseddata)
		
		if(normalisation=="mas5")exprdata<-log(exprdata,2)
		annotation<-.convertEnsembl(platform,values=rownames(exprdata))
		datagenes<-.combineEnsembl(exprdata,annotation,genelist)
		d2<-.fitlms(normalised=datagenes,customfactors=customfile,experiment=experimentname,case=case,statistic=statistic)
	}
	print('Fitted linear models: ')
	print(colnames(d2$ranklist))

	return(d2)
}
