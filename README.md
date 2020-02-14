Integrating GEO Transcriptone Data in Random-combining to Find Differentiated Genes in R (ITRDR)
==========
ITRDR是Ｒ平台下的基于不同来源的转录组大数据的整合工具，适用于ncbi geo等转录组大数据的整合及其差异基因筛选和特定通路基因筛选等工作，
筛选到的基因可直接用于实验验证准确率较高。

## Authors

[刘文粟] (Liu Wensu)

[周晓北] (Zhou Xiaobei)



## Usage
读取GEO数据并形式上合并
```r
library(data.table)
dir <- "~/ITRDR/data/sample_data"
setwd(dir)
f <- list.files(pattern="^GSE",full.names=TRUE)
f1 <- lapply(f,fread)
nms <- gsub(".csv","",basename(f))
f1 <- mapply(function(x,y) {nm <- names(x)
                            idx <- !nm=="gene"
                            names(x)[idx] <- paste0(y,"_",nm[idx])
                            x},x=f1,y=nms,SIMPLIFY=FALSE)

f1 <- Reduce(function(...) merge.data.table(...,all=TRUE),f1)
#f1 <- f1[grep("[0-9][0-9]-[A-Z][a-z][a-z]|^__",gene,invert=TRUE)]
f1[is.na(f1)] <- 0
fwrite(f1,"combine.csv")
```

读入文献中有报道的标志性基因和已知通路
```r
marker_gene<-c("AR","CD46","STAT5A","STAT5B","NOTCH1","FOXA1",
   "KDM6A","ASXL1","KMT2A","KMT2D","WDR5","ASH2L","EZH2",
   "KDM1A","NR3C1","PIK3CB","AKT1","MTOR","KDM5D","C17ORF49")
power_gene<-sqrt(c(16,8,5,5,8,6,6,6,10,10,10,10,8,8,12,8,8,8,8,8)) 
target_GO<-c("GO:0000380","GO:0016310","GO:0030521","GO:0008083","GO:0006914")
library("org.Hs.eg.db")
go <- as.data.table(org.Hs.egGO)
geneinfo <- as.data.table(org.Hs.egSYMBOL)
go <- merge(go,geneinfo,all.x=TRUE)
gene_go <- go[go_id%in%target_GO,unique(symbol)]
```
寻找研究对象的markers
```r
topic <- "and prostate cancer"
library(data.table)
x <- fread("combine.csv")
gene <- x[,gene]
x[,gene:=NULL]
y <- list()
for(i in seq(gene))
{
   tmp <- getCount(gene[i])
   if(tmp>0)
   {
      xi <- paste(gene[i],topic)
      y[[i]] <- getCount(xi)
   }
   else
   {
      y[[i]] <- 0
   } 
   cat("iteration = ",i,"\n")
}
y <- unlist(y)
res <- data.table(gene=gene,hits=y)
fwrite(res,"hit.csv")
```