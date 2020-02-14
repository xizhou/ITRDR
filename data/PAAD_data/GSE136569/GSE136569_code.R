library(data.table)
f <- list.files(pattern="^GSM")
mm <- rep("",length(f))
mm[grep("N",f)] <- paste("Normal_",seq(grep("N",f)))
mm[grep("T",f)] <- paste("Tumor_",seq(grep("T",f)))
f1 <- lapply(f,fread)
names(f1) <- mm
f1 <- lapply(f1,function(x) x[,.(RPB,Symbol)])
f1 <- lapply(f1,function(x) x[,RPB:=sum(RPB),by=Symbol])
f1 <- lapply(f1,unique)
f1 <- lapply(f1,function(x) x[!Symbol=="",])
f1 <- mapply(function(x,y) {x <- x[,(y):=RPB]
                       x[,RPB:=NULL]
                       x},x=f1,y=mm,SIMPLIFY=FALSE)
tpm <- f1[[1]]
for(i in 2:length(f1))
   tpm <- merge(tpm,f1[[i]],all=TRUE)
tpm <- tpm[grep(":|-",Symbol,invert=TRUE),]
tpm[is.na(tpm)] <- 0
tpm[,gene:=Symbol]
tpm[,Symbol:=NULL]
gn <- tpm[,gene]
tpm[,gene:=NULL]

gtf <- fread("Homo_sapiens.GRCh38.99.gtf",header=FALSE)
l <- gtf[,V5-V4]
g <- gtf[,V9]
g <- strsplit(g,";")
g <- sapply(g,function(x) x[grep(" gene_name",x)])
g <- gsub(" gene_name \"|\"","",g)
ano <- data.table(length=l,gene=g)
ano <- ano[,max(length),by=gene]
l <- ano[,V1][match(gn,ano[,gene])]

library(tpmToFpkm)

fpkm <- tpmToFpkm(tpm[!is.na(l),],l[!is.na(l)])
fpkm <- as.data.table(fpkm)
fpkm[,gene:=gn[!is.na(l)]]

out <- fpkm
out[,meanN:=rowMeans(.SD+1e-8),.SDcols=grep("Normal",names(out),value=T)]
out[,meanT:=rowMeans(.SD+1e-8),.SDcols=grep("Tumor",names(out),value=T)]
out[,fold_T_NM:=(meanT/meanN)]
fwrite(out,"GSE136569.csv")