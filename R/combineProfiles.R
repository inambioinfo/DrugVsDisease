combineProfiles<-function(data){
	if(!is.list(data)){stop('Please provide list of R objects for data')}
	refdata<-data[[1]]
	norm<-list(length=length(data))
	for(i in 2:length(data)){
		#pairwise do the comparison to the ref data
		out<-.mrs(refdata,data[[i]])	
		norm[[i]]<-out$y
	}	
	norm[[1]]<-out$x
	final<-norm[[1]]
	for(i in 1:length(data)){
		temp<-norm[[i]]
		final<-cbind(final,temp[rownames(norm[[1]]),])
	}
	return(final)
}
