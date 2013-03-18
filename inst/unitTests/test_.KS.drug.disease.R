test_.KS.drug.disease<-function(){
	library(RUnit)
	library(DrugVsDiseasedata)
	data(diseaseRL,package="DrugVsDiseasedata")
	d<-diseaseRL[,1:3]
	checkTrue(sum(diag(DrugVsDisease:::.KS.drug.disease(d,d)))==3)
	checkTrue(all(DrugVsDisease:::.KS.drug.disease(d,d)<=1))
	checkTrue(all(DrugVsDisease:::.KS.drug.disease(d,d)>=-1.03))
	
}
