test_.getESscores<-function(){
	library(RUnit)
	library(DrugVsDiseasedata)
	data(diseaseRL,package="DrugVsDiseasedata")
	d<-diseaseRL[,1]
	d<-as.matrix(d)
	d1<-DrugVsDisease:::.ranknamemat(d,2,TRUE)
	d2<-DrugVsDisease:::.ranknamemat(d,2,FALSE)
	checkTrue(!is.null(colnames(DrugVsDisease:::.getESscores(d1,d2,drugPRLs=diseaseRL[,1:20]))))
	
}
