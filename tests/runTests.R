require("DrugVsDisease")||stop("unable to load DrugVsDisease")
library(BiocGenerics)
library(RUnit)
BiocGenerics:::testPackage("DrugVsDisease")