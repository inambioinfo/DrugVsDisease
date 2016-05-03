selectrankedlists <-
function(ranklist,colsinc){
	#get the selected ranklists/pvalues as chosen in colsinc
	seldata<-vector("list",length=2)
	cnames<-colnames(ranklist[[1]])

		seldata[[1]]<-as.matrix(ranklist$ranklist[,colsinc])
		seldata[[2]]<-as.matrix(ranklist$pvalues[,colsinc])
		colnames(seldata[[1]])<-cnames[colsinc]
    colnames(seldata[[2]])<-cnames[colsinc]
	
	names(seldata)<-names(ranklist)
	return(seldata)
	
}
