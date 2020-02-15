Integrating GEO Transcriptone Data in Random-combining to Find Differentiated Genes in R (ITRDR)
==========
ITRDR是Ｒ平台下的基于不同来源的转录组大数据的整合工具，适用于ncbi geo等转录组大数据的整合及其差异基因筛选和特定通路基因筛选等工作，
筛选到的基因可直接用于实验验证准确率较高。

## Authors

[刘文粟] (Liu Wensu)

[周晓北] (Zhou Xiaobei)

## <a name="install"></a> Prerequisites and installation

### <a name="prerequisites"></a> Prerequisites
- **R (>=3.5.3)** 

### <a name="installlation"></a> Installlation
Temporally, use:
```r
source("~/ITRDR/R/ITRDR.R")
```

## Usage
### <a name="process"></a> Process data
The input data should be formated as following: 

    gene	pval_1
    NECAP2	0.358843944
    CLIC4	0.030455456
    DBT	0.292270083
    C1orf21	0.000382696
    COP1	0.020345337
    PRUNE1	0.325788738
    LIN9	0.017277068

The first column is gene symbol (with column name's "gene"), the other columns should contain the result of selected group comparison (pvalue (named as pval\_x) or 
fold change (fc\_x)).

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
The result should be like:

    gene	GSE55030_p_val	GSE71797_fc_1	GSE71797_fc_2	GSE82223_fc_1	GSE86532_fc_1	GSE92574_fc_1	GSE110903_fc_1
    ISG15	0.029289252	0.475583449	3.077468801	0	0.126152927	0	0.484946133
    AGRN	0.020108192	9.767047772	0.460773615	0	0	0	1.677513801
    RP11-465B22.3	0	0	0	0	0	0	0.125557896
    RP11-465B22.8	0	0	0	0	0	0	0.395281318
    B3GALT6	0.391504072	2.840062329	0.312959305	0	0	0.334790739	1.010640313

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