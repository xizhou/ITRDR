library(stringr)
library(biomaRt)
library(affy)
library(data.table)
#f <- list.files(pattern="CEL")
f <- fread("sample.txt")
a=ReadAffy(filenames=f[,sample])
sampleNames(a) <- f[,treatment]
x <- expresso(a,bg.correct=TRUE,bgcorrect.method="mas",normalize.method="constant",
   pmcorrect.method="pmonly",summary.method="medianpolish")

x <- exprs(x)
rns <- rownames(x)
x <- as.data.table(x)
#x[,gene:=str_sub(rns,2,-1)]
x[,gene:=rns]
mart <- useMart("ensembl", "hsapiens_gene_ensembl")   #use human data
filters <- listFilters(mart)
y <- getBM(attributes=c("hgnc_symbol","affy_hugene_1_0_st_v1"),mart=mart)   ##change term according to affy type
y <- as.data.table(y)
y[,gene:=affy_hugene_1_0_st_v1]
y[,affy_hugene_1_0_st_v1:=NULL]
y[,gene:=as.character(gene)]
out <- merge(x,y,all.x=TRUE)
out1 <- split(out,by="gene")
out1 <- lapply(out1,"[",1L)
out <- rbindlist(out1)
out[,gene:=NULL]
out[,gene:=hgnc_symbol]
out[,hgnc_symbol:=NULL]
out <- out[(!is.na(gene))&(!gene==""),]
cols <- names(out)[!names(out)=="gene"]
out[,(cols):=2^(.SD),.SDcols=cols]
out <- out[,lapply(.SD,mean),by=gene]
out[,meanN:=rowMeans(.SD+1e-8),.SDcols=grep("Normal",names(out),value=T)]
out[,meanT:=rowMeans(.SD+1e-8),.SDcols=grep("Tumor",names(out),value=T)]
out[,fold_T_NM:=(meanT/meanN)]
fwrite(out,"GSE62452.csv")

