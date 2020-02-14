Integrating GEO Transcriptone Data in Random-combining to Find Differentiated Genes in R (ITRDR)
==========
ITRDR是Ｒ平台下的基于不同来源的转录组大数据的整合工具，适用于ncbi geo等转录组大数据的整合及其差异基因筛选和特定通路基因筛选等工作，
筛选到的基因可直接用于实验验证准确率较高。

## Authors

[刘文粟] (Liu wensu)

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