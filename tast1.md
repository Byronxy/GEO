---
title: "prostate cancer"
author: "Yi Xiong"
date: "2018年10月16日"
output: html_document
---
This is homework for biotree
### To reproduce the result from https://www.jci.org/articles/view/96060/figure/1 c-f
TRAF4 gene
##Figure1f

```{r}
setwd("E:/bioinformatics/GEO/prostate cancer")
library(GEOquery)
```
获取GEO数据
```{r}
#getGEO函数是关键，通常过去GSE number
#gset <- getGEO('GSE3325',
#                 AnnotGPL = F,
#                 getGPL = T)
#save(gset,file = 'GSE3325.gset.Rdata')
load(file = 'GSE3325.gset.Rdata')
```
获取基因表达矩阵，样本矩阵。
```{r}
if(F){
gset = gset[[1]]
expreSet = exprs(gset)
pdata = pData(gset)
#grouplist = as.character(pdata$description)
dim(expreSet) #54675    19
#根据描述取出 normal,primary,metastatic cancer
index = c(1:4,7:11,14:17)
expreSet = expreSet[,index]
expreSet[1:5,]
grouplist = c(rep("normal",4),rep("primary",5),rep("metastatic",4))
table(grouplist)
}
```
接下来处理探针
```{r}
if(F){
GPL = gset@featureData@data #从getGEO中下载了平台，获得基因和探针之间的关系
#其中ID为探针，Gene Symbol为基因名
ids = GPL[,c(1,11)]
#去除掉没有gene symbal的探针
ids = ids[ids$`Gene Symbol`!="",]
table(sort(table(ids$`Gene Symbol`)))
tail(sort(table(ids$`Gene Symbol`)))
#一个探针对应了多个基因，干脆不处理算了，或者把第二个基因去掉
a<-strsplit(as.character(ids[,2]), " /// ")
tmp <- mapply(cbind, ids[,1], a ) 
#这个办法最好
tmp = plyr::ldply(tmp, rbind)
tmp[1:5,]
rownames(tmp) = make.unique(tmp[,1])
tmp = tmp[,-1]
colnames(tmp)=c("ID","Symbal")
ids = tmp
#这一步根据ids去除没有注释的探针
expreSet = expreSet[ rownames(expreSet) %in% ids[ , 1 ], ]
ids = ids[ match(rownames(expreSet), ids[ , 1 ] ), ]
dim( expreSet )
dim( ids )
#多个探针对应一个基因的取最大值
tmp = by(expreSet,ids$Symbal,function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
expreSet=expreSet[rownames(expreSet) %in% probes ,]
rownames(expreSet ) = ids[match(rownames(expreSet), ids[,1]), 2 ]
#获取表达量
expreSet = log(expreSet)
}

#save(expreSet,grouplist,file = 'GSE3325.expreSet.Rdata')
load(file = 'GSE3325.expreSet.Rdata')
#plot
TRAF4 = expreSet[rownames(expreSet)=="TRAF4",]
plotdata = data.frame(grouplist = grouplist, values = TRAF4, stringsAsFactors = F)
plotdata$grouplist = factor(plotdata$grouplis, levels = c("normal", "primary", "metastatic"))
library(ggstatsplot)
p1 = ggbetweenstats(data = plotdata, x = grouplist,  y = values)
library(ggpubr)
my_comparisons <- list(c("normal", "primary"), c("primary", "metastatic"), c("normal", "metastatic"))
ggboxplot(plotdata, x ="grouplist", y = "values",
          color = "grouplist", palette = "jco",
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons)
```
Additional analysis
```{r}
library("factoextra")
df = scale(t(expreSet))
rownames(df) = make.unique(grouplist)
dist <- dist(df)
hc = hclust(d = dist, method = "average")
fviz_dend(hc, cex = 0.5)
grp <- cutree(hc, k = 3)
head(grp, n = 3)
table(grp)
fviz_dend(hc, k = 3, # Cut in four groups
cex = 0.5, # label size
k_colors = c("#2E9FDF", "#00AFBB", "#E7B800"),
color_labels_by_k = TRUE, # color labels by groups
rect = TRUE # Add rectangle around groups
)
```

Figure1c chandran and yu
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6919
GSE6919
```{r}
# gset <- getGEO('GSE6919',
#                 AnnotGPL = F,
#                 getGPL = T)
# save(gset,file = 'GSE6919.gset.Rdata')
load(file = 'GSE6919.gset.Rdata')
```
```{r}
gset = gset[[1]]
expreSet = exprs(gset)
pdata = pData(gset)
dim(expreSet)
#phenotype data
grouplist = as.character(pdata$description)
table(grouplist)
```
```{r}
if(F){
  GPL = gset@featureData@data #从getGEO中下载了平台，获得基因和探针之间的关系
#其中ID为探针，Gene Symbol为基因名
ids = GPL[,c(1,11)]
#去除掉没有gene symbal的探针
ids = ids[ids$`Gene Symbol`!="",]
table(sort(table(ids$`Gene Symbol`)))
tail(sort(table(ids$`Gene Symbol`)))
#一个探针对应了多个基因，干脆不处理算了，或者把第二个基因去掉
a<-strsplit(as.character(ids[,2]), " /// ")
tmp <- mapply(cbind, ids[,1], a ) 
#这个办法最好
tmp = plyr::ldply(tmp, rbind)
tmp[1:5,]
tmp = na.omit(tmp)
rownames(tmp) = make.unique(tmp[,1])
tmp = tmp[,-1]
colnames(tmp)=c("ID","Symbal")
ids = tmp
#这一步根据ids去除没有注释的探针
expreSet = expreSet[rownames(expreSet) %in% ids[, 1], ]
ids = ids[match(rownames(expreSet), ids[ , 1]), ]
dim( expreSet )
dim( ids )
#多个探针对应一个基因的取最大值
tmp = by(expreSet,ids$Symbal,function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character(tmp)
expreSet=expreSet[rownames(expreSet) %in% probes ,]
rownames(expreSet ) = ids[match(rownames(expreSet), ids[,1]), 2 ]
#获取表达量
expreSet = log(expreSet)
}

#save(expreSet,grouplist,file = 'GSE6919.expreSet.Rdata')
load(file = 'GSE6919.expreSet.Rdata')
```
```{r}
#plot
TRAF4 = expreSet[rownames(expreSet)=="TRAF4",]
plotdata = data.frame(grouplist = grouplist, values = TRAF4, stringsAsFactors = F)
library(ggstatsplot)
p1 = ggbetweenstats(data = plotdata, x = grouplist,  y = values)
#似乎没有这个基因
```

```{r}
# gset <- getGEO('GSE6919',
#                 AnnotGPL = F,
#                 getGPL = T)
# save(gset,file = 'GSE6919.gset.Rdata')
load(file = 'GSE6919.gset.Rdata')
```

