test_.getESscores<-function(){
	library(RUnit)
	library(DvDdata)
	data(drugRL,package="DvDdata")
	d<-drugRL[,1]
	d<-as.matrix(d)
	d1<-DrugVsDisease:::.ranknamemat(d,2,TRUE)
	d2<-DrugVsDisease:::.ranknamemat(d,2,FALSE)
	checkTrue(!is.null(colnames(DrugVsDisease:::.getESscores(d1,d2,drugPRLs=drugRL[,1:20]))))
	
}