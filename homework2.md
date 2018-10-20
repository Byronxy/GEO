---
title: "STAT115 Homework 2 microarray analysis"
author: "(Yi Xiong)"
date: "`r Sys.Date()`"
---

# Part II: Microarray Clustering and Classification

The sample data is in file "taylor2010_data.txt" included in this
homework. This dataset has expression profiles of 23974 genes in 27
normal samples, 129 primary cancer samples, 18 metastasized cancer
samples, and 5 unknown samples. Assume the data has been normalized and
summarized to expression index.

```{r loadtaylor}
setwd("E:/bioinformatics/lecture/homework2")
taylor <- as.matrix(read.csv("taylor2010_data.txt", sep="\t",row.names=1))
index_normal <- grepl("N.P", colnames(taylor))
index_primary <- grepl("P.P", colnames(taylor))
index_met <- grepl("M.P", colnames(taylor))
n_normal <- sum(index_normal);
n_primary = sum(index_primary);
n_met = sum(index_met);

# class label (design vector)
taylor_classes = c(rep(0,n_normal), rep(1,n_primary), rep(2,n_met));

# train (known type samples), and test (unknown type samples)
train <- taylor[,1:174];
test <- taylor[,175:179];

# colors for plotting
cols = c(taylor_classes+2, 1,1,1,1,1)

tumortype_class <- factor(taylor_classes, levels = 0:2,
                          labels = c("Normal", "Primary", "Metastasized"))

train_samps <- 1:174
test_samps <- 175:179
```

**1. For the 174 samples with known type (normal, primary, metastasized), use LIMMA to find the differentially expressed genes with fold change threshold 1.5, and adjusted p-value threshold 0.05. How many differentially expressed genes are there?**
Hint: the design vector `taylor_classes` consists of type indicator for the 174 samples. For example, 0 for normal, 1 for primary, and 2 for metastasized.

```{r limma3}
#differentially expressed genes
library(limma)
#design matrix
design.mat <- model.matrix(~0+factor(tumortype_class))
colnames(design.mat) <- c("Normal", "Primary", "Metastasized")
head(design.mat)
# contrast matrix
contrast.mat <- makeContrasts(T1 = Primary-Normal, T2 = Metastasized-Normal, T3 = Metastasized-Primary, levels = design.mat)
contrast.mat
fit <- lmFit(train, design.mat)
fit2 <- contrasts.fit(fit, contrast.mat)
fit3 <- eBayes(fit2)
summary(decideTests(fit3, adjust.method = "fdr", p.value = 0.05, lfc = log2(1.5)))
result <- decideTests(fit3, adjust.method = "fdr", p.value = 0.05, lfc = log2(1.5))
taylor.expressed <- taylor[result[,1]!=0|result[,2]!=0|result[,3]!=0,]#select all differentially expressed genes |= or
taylor.expressed[1:20,1:5]
dim(taylor.expressed) #there are 3822 differentially expressed genes
#把3个组的差异表达基因分别提取出来
PvsN <- topTable(fit3, coef = 1, number = Inf, p.value = 0.05, lfc = log2(1.5))
MvsN <- topTable(fit3, coef = 2, number = Inf, p.value = 0.05, lfc = log2(1.5))
MvsP <- topTable(fit3, coef = 3, number = Inf, p.value = 0.05, lfc = log2(1.5))
length(rownames(PvsN))
length(rownames(MvsN))
length(rownames(MvsP))
```
画一个韦恩图
```{r}
library(VennDiagram)
library(RColorBrewer)
venn <- venn.diagram(list(PvsN = rownames(PvsN), MvsN=rownames(MvsN),MvsP =rownames(MvsP))
             ,filename = "venn_three.tiff",#filename 
             fill=c(brewer.pal(3, "Set1")),
             col = "white")#three elements
```
![Caption for the picture.](venn_three.tiff)
批次处理。做一个聚类分析，感觉是有批次效应的。但是没有给出批次的信息，所以暂时不做校正
```{r}
t.taylor.expressed <- t(taylor.expressed) #转置
group <- c(rep("Normal",n_normal),rep("Primary",n_primary),
           rep("Metastasis",n_met),rep("Unknow",5))
rownames(t.taylor.expressed) <- make.unique(group)
dim(t.taylor.expressed)
deg.clusters <- hclust(dist(t.taylor.expressed), method = 'average')
dend <- as.dendrogram(deg.clusters)
library(dendextend)
colorCodes <- c(Normal = "red", Primary = "green", Metastasis = "blue", Unknow = "black")
labels_colors(dend) <- colorCodes[group][order.dendrogram(dend)]
plot(dend, cex = 0.8, main = "clustering")
```

**2. Draw k-means clustering on all samples using the differentially expressed genes. Do the samples cluster according to disease status?**

```{r kmeans}
set12 <-union(rownames(PvsN),rownames(MvsN))
set123 <-union(set12,rownames(MvsP))
length(set123)
library(sva)
edata <- taylor[set123,]
kClust <- kmeans(t(edata), centers=3, nstart = 25, iter.max = 20)
kClusters <- as.factor(kClust$cluster)
taylor_classes = c(rep("Normal",n_normal), 
                   rep("Primary",n_primary), 
                   rep("Metastasized",n_met),
                   rep("Unknown",5))
taylor_classes <- factor(taylor_classes,
                         levels = c("Normal", "Primary", "Metastasized","Unknown"))
plotdata <- data.frame(taylor_classes,kClusters)
library(ggplot2)
ggplot(plotdata, aes(x = taylor_classes, y = kClusters, colour = kClusters)) +
  # 散点图函数
  geom_point(position="jitter")
```
我们发现，癌和正常组有一部分是重叠的，而未知样品有三个属于转移样品，剩余的两个要么属于正常，要么就是肿瘤。 这里分不开是有原因的，我们从韦恩图上看出，癌和癌旁的从有一部分差异基因就是重合的。
突发奇想，我们把三个角上的基因拿出来，这样理论上对三个群体都会有很好的分类。
提取数据
```{r}
PvsN_unique <- setdiff(rownames(PvsN),union(rownames(MvsP),rownames(MvsN)))
MvsN_unique <- setdiff(rownames(MvsN),union(rownames(MvsP),rownames(PvsN)))
MvsP_unique <- setdiff(rownames(MvsP),union(rownames(PvsN),rownames(MvsN)))
all_unique <-union(union(PvsN_unique,MvsN_unique),MvsP_unique)
```
绘图
```{r}
library(sva)
edata <- taylor[all_unique,]
kClust <- kmeans(t(edata), centers=3, nstart = 25, iter.max = 20)
kClusters <- as.factor(kClust$cluster)
taylor_classes = c(rep("Normal",n_normal), 
                   rep("Primary",n_primary), 
                   rep("Metastasized",n_met),
                   rep("Unknown",5))
taylor_classes <- factor(taylor_classes,
                         levels = c("Normal", "Primary", "Metastasized","Unknown"))
plotdata <- data.frame(taylor_classes,kClusters)
library(ggplot2)
ggplot(plotdata, aes(x = taylor_classes, y = kClusters, colour = kClusters)) +
  geom_point(position="jitter")
```
这时候我们发现，没能解决问题，还带来新的问题，转移的和肿瘤分不清楚了，原因是我们去掉的那1119个基因是专门用来区分肿瘤和转移的。

**3. Do PCA biplot on the samples with differentially expressed genes genes, and use 4 different colors to distinguish the 4 types of samples (normal, primary, metastasized and unknown). Do the samples from different groups look separable?**
Hint: use the PCA ggplot R code, also function legend is useful. (http://docs.ggplot2.org/0.9.3.1/geom_point.html)

```{r pca-biplot}
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
color <- c(brewer.pal(4,"Set1"))
pc.all = prcomp(t.taylor.expressed, scale. = T)
sumpca <- summary(pc.all)
xlab <- paste("PC1","(",round(summary(pc.all)$importance[2,1],2),"% explained variance)",sep = "")
ylab <- paste("PC2","(",round(summary(pc.all)$importance[2,2],2),"% explained variance)",sep = "")
group <- as.factor(group)
unknow <- rep("Unknow",5)
data_x <- data.frame(pc.all$x[,1:2]) #extract PC1 and PC2
rownames(data_x) <- make.unique(c(rep("Normal",n_normal),rep("Primary",n_primary),
                      rep("Metastasis",n_met),rep("Unknow",5)))
ggplot(data_x,aes(PC1,PC2))+geom_point(aes(color=group),size=4,alpha=0.6)+
  scale_color_manual(values = color)+geom_text_repel(data = data_x[175:179,],aes(label = unknow),
                                                     box.padding=unit(0.5, "lines"))+
  stat_ellipse(data = data_x[1:174,], aes(PC1,PC2, fill = group[1:174]), geom = "polygon", alpha = 1/2, level = 0.95, show.legend = F)+
  xlab(xlab)+ylab(ylab)
```

**4. FOR GRADUATES: What percent of variation in the data is captured in the first two principle components? How many principle components do we need to capture 85% of the variation in the data?**
R Hint: use function prcomp.

```{r percent-variance}
std_dev <- pc.all$sdev #compute standard deviation of each principal component
pr_var <- std_dev^2 #compute variance
pr_var[1:10] #check variance of first 10 components
prop_varex <- pr_var/sum(pr_var) #proportion of variance explained
library(ggplot2) #PC解释方差的累积曲线
tmp <- data.frame(n <- c(1:179), cumsum <- cumsum(prop_varex))
tmp[round(tmp$cumsum....cumsum.prop_varex.,2)==0.85,]
ggplot(tmp)+geom_point(aes(x = n, y = cumsum), 
                       xlab = "Principal Component",
                       ylab = "Cumulative Proportion of Variance Explained")+geom_hline(yintercept=0.85, linetype="dashed", color = "red")
```

**5. Based on the PCA biplot, can you classify the 5 unknown samples?  Put the PCA biplot in your HW write-up, and indicate which unknown sample should be classified into which known type (normal, primary, metastasized). Do you have different confidence for each unknown sample?**
3个样品肯定是转移的，另外两个可能是正常或者原发

**6. FOR GRADUATES: Use PCA on all samples and all the genes (instead of the differentially expressed genes) for sample classification. Compare to your results in the previous question. Which PCA plot looks better and why?**

```{r pca-all}
t.taylor <- t(taylor)
pca = prcomp(t.taylor, scale. = T)
library(ggbiplot)
ggbiplot(pca, choices = 1:2,groups = group,
                  ellipse.prob = 0.5, circle = F, obs.scale = 1, var.scale = 1,
                  var.axes = F, ellipse = T)+theme_bw()
```
肯定有人会不理解，为什么要先算差异基因呢，最起码有两个原因： 第一，也是最重要的，差异基因在本质上校正了批次效应，排除了很多背景。 第二，是我自己发现的，实际上不重要，SVM还有回归这些算法不支持大量基因。 感兴趣的翻看到公众号的SVM那个帖子看看就知道了，他也是要接受差异基因的。

**7. Run KNN (try K = 1, 3 and 5) on the differential genes and all the samples, and predict the unknown samples based on the 174 labeled samples. Hint: use the library `class` and function `knn`.**

```{r knn}
library(class)
train <- taylor[set123,1:174];
test <- taylor[set123,175:179];
knn(t(train),t(test),tumortype_class,k=1,prob = T)
knn(t(train),t(test),tumortype_class,k=3,prob = T)
knn(t(train),t(test),tumortype_class,k=5,prob = T)
```

**8. Run SVM (try a linear kernel) on the differential genes and all the samples, and predict the unknown samples based on the 174 labeled samples. Hint: use the library `e1071` and function `svm`.**

```{r svm}
library(e1071)
edata <- taylor[set123,]
svm_data <- t(edata)
taylor_classes = c(rep(0,n_normal), rep(1,n_primary), rep(2,n_met))
tumortype_class <- factor(taylor_classes, levels = 0:2,
                          labels = c("Normal", "Primary", "Metastasized"))
m <- svm(svm_data[train_samps,], tumortype_class,kernel = "linear")
m
predict(m,svm_data[test_samps,])

```
**9. FOR GRADUATES: Implement a 3-fold cross validation on your SVM classifier, based on the 174 samples with known labels. What is your average (of 3) classification error rate on the training data?**
使用K-fold交叉验证
（将原始数据分成 K 组 (一般是均分), 将每个子集数据分别做一次验证集, 其余的 K-1 组子集数据作为训练集, 这样会得到 K 个模型, 用这 K 个模型最终的验证集的分类准确率的平均数作为此 K-CV 下分类器的性能指标. K 一般大于等于 2, 实际操作时一般从 3 开始取。）
```{r}
train_data = as.data.frame(svm_data[train_samps,])
train_data$tumortype_class = tumortype_class
svm_tune <- tune(svm, tumortype_class ~ ., data = train_data,
                 kernel = "linear",
                 ranges = list(cost = c(0.01, 0.1, 1, 10)),
                 tunecontrol = tune.control(cross = 3))
svm_tune
```
模型构建和预测,注意这里预测的是训练数据
```{r}
m <- svm(svm_data[train_samps,], tumortype_class,kernel = "linear",cost = 0.01)
svm.pred <- predict(m,svm_data[train_samps,])
```
计算错误率
```{r}
tab <- table(pred = svm.pred, true = tumortype_class)
tab
classification_error <- 1- sum(svm.pred == tumortype_class)/length(svm.pred)
classification_error
```
发现是0，也就是说这个模型拟合的很好。
在这个基础上继续预测
```{r}
predict(m,svm_data[test_samps,])
```

总结:Kmeans和PCA靠看，SVM和KNN靠算，靠看容易看花眼，靠算长留天地间。

但是，还是有问题，在这里SVM把自己拟合的十分美好，存在过度拟合的嫌疑，这样就失去了普适性。

不过，不要担心，lasso回归可以解决这个问题，如果他自己拟合的太好，他就会惩罚自己。 我们还用过lasso回归来做生存模型构建呢。
