test_.KS.drug.disease<-function(){
	library(RUnit)
	library(DvDdata)
	data(drugRL,package="DvDdata")
	d<-drugRL[,1:3]
	checkTrue(sum(diag(DrugVsDisease:::.KS.drug.disease(d,d)))==3)
	checkTrue(all(DrugVsDisease:::.KS.drug.disease(d,d)<=1))
	checkTrue(all(DrugVsDisease:::.KS.drug.disease(d,d)>=-1))
	
}