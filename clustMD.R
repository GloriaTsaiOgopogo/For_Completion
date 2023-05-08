file <- file.choose()
data <-read.table(file, header = T, sep = "\t")

head(data)
dim(data)
## [1] 982  48
data[complete.cases(data),]
dim(data[complete.cases(data),])
## [1] 19 48


d1 <-subset(data, select = -c(c(门诊号,姓名, 手机号)))
d1 <-subset(data, select = -c(c(ID)))

n <- 1
x <-vector()

for(i in 1:ncol(d1)) {
  x[n] <- sum(complete.cases(d1[,i]))
  n <- n + 1
}
names(x) <- colnames(d1)

sort(x)

d1 <- subset(d1, select = -c(c(F, LLDL, HbA1C)))

dim(d1[complete.cases(d1), ])
# [1] 38 38


d1 <- subset(d1, select = -c(E2))

dim(d1[complete.cases(d1), ])
complete.d1 <- d1[complete.cases(d1), ]

str(data)

data$B超 <-as.factor(data$B超)
data$月经模式 <- as.factor(data$月经模式)
data$家族史 <- as.factor(data$家族史)
str(data)


# 调和曲线
unison<-function(x){
   # x is a matrix or data frame of data
   if (is.data.frame(x)==TRUE)
      x<-as.matrix(x)
   t<-seq(-pi, pi, pi/30)
   m<-nrow(x); n<-ncol(x)
   f<-array(0, c(m,length(t)))
   for(i in 1:m){
      f[i,]<-x[i,1]/sqrt(2)
      for( j in 2:n){
          if (j%%2==0) 
             f[i,]<-f[i,]+x[i,j]*sin(j/2*t)
          else
             f[i,]<-f[i,]+x[i,j]*cos(j%/%2*t)
      } 
   }
   plot(c(-pi,pi), c(min(f),max(f)), type="n", 
        main="The Unison graph of Data",
        xlab="t", ylab="f(t)")
  for(i in 1:m) lines(t, f[i,] , col=i)
}

unison(data)

library(lattice)
parallel(~data[1:37], data, 
    horizontal.axis = FALSE, scales = list(x = list(rot = 90)))

data <-subset(data, select = -c(c(FSH,DHEA.S,LH,PRL,F)))

library(cluster)
kc=pam(data,5)
#kc=pam(data,3,cluster.only=TRUE)
print(kc)
summary(kc)
plot(kc)

library(factoextra)
df<-scale(data)
fviz_nbclust(df, kmeans, method = "wss") + geom_vline(xintercept = 4, linetype = 2)

set.seed(12345)
km_result <- kmeans(df, 4, nstart = 24)
print(km_result)
dd <- cbind(df, cluster = km_result$cluster)
head(dd)
dd <- as.data.frame(dd)
table(dd$cluster)

fviz_cluster(km_result, data = df,
             palette = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),
             ellipse.type = "euclid",
             star.plot = TRUE, 
             repel = TRUE,
             ggtheme = theme_minimal()
)

result <- dist(df, method = "euclidean")

##############################################
library(clustMD)

data <-read.table(file.choose(), header = T, sep = "\t")

#data$B超 <-as.numeric(data$B超) 
#data$月经模式 <- as.numeric(data$月经模式)
#data$家族史 <- as.numeric(data$家族史)
str(data)

## Order variables (Continuous, ordinal, nominal)
data <- data[, c(1:4,12:40,7:10,5,11,6)]
str(data)

#data <- as.matrix(data)

# Standardise continuous variables
data[, c(1:33)] <- scale(data[, c(1:33)])


# Start categorical variables at 1 rather than 0
data[, 34:40] <- data[, 34:40] + 1

head(data)

res1 <- clustMD(X = data, G = 2, CnsIndx = 33, OrdIndx = 40, Nnorms = 100000,
    MaxIter = 500, model = "EVI", store.params = FALSE, scale = TRUE, 
    startCL = "kmeans", autoStop= TRUE, ma.band=30, stop.tol=0.0001)

res2 <- clustMD(X = data, G = 2, CnsIndx = 33, OrdIndx = 40, Nnorms = 1000000,
    MaxIter = 500, model = "VII", store.params = FALSE, scale = TRUE, 
    startCL = "hc_mclust", autoStop= TRUE, ma.band=30, stop.tol=0.0001)

library(snow)
res <- clustMDparallel(X = data, G = 1:2, CnsIndx = 33, OrdIndx = 40, Nnorms = 200000,
	MaxIter = 500, models = c("EVI", "EII", "VII"), store.params = FALSE, scale = TRUE,
	startCL = "hc_mclust", autoStop= TRUE, ma.band=30, stop.tol=0.0001)

abs(res1$cl - res2$cl) == 1
table(abs(res1$cl - res2$cl) == 1)
FALSE  TRUE 
   49    76 
   
cluster_tag <- res2$cl

res1$BIChat
res2$BIChat

data <-read.table(file.choose(), header = T, sep = "\t")

data$cl <- cluster_tag 

write.table(data,file.choose(), sep = "\t", row.names = F)


## 选取BIC小的数据

## CnsIndx 连续变量的个数
## OrdIndx 连续变量、二进制变量和有序变量的个数（不含名义变量）

summary(res)
res$cl  ### 预测的分类
plot(res)
===========================================================================

bi <- subset(data, select = c(c(家族史,月经模式,B超)))
nu <- subset(data, select = -c(c(家族史,月经模式,B超)))

tmp <- c() 
for(i in 1:ncol(nu)){
  # 使用shapiro.test函数检验每个变量是否与正太分布具有显著性差异
  p_value <- apply(nu,2,shapiro.test)[[i]]$p.value
  # 输出结果变量名和p值。
  tmp[i] <- p_value
}
names(tmp) <- colnames(nu)
sort(tmp)
sort(tmp > 0.01)

library(reshape2)
library(ggplot2)
library(ggpubr)

d <-melt(data,id.var="cl")

## p<-ggplot(data = d, aes(x = factor(cl),y = value))+geom_boxplot(aes(fill=variable))

ggplot(data = data, aes(x=oocyteType,y=number))+geom_boxplot(aes(fill=fociType))

for i in names(data){
	s1 <- ggplot(data = d[d$variable == i,], aes(x = factor(cl), y = value, fill = variable)) + geom_boxplot() + stat_compare_means(method = "t.test")
}

s1 <- ggplot(data = d[d$variable == names(data)[40],], aes(x = factor(cl), y = value, fill = variable))
s1 <- s1 + geom_boxplot(fill = c('#FFCC00','#FF9900')) + geom_point(position = position_jitterdodge()) + stat_compare_means(method = "t.test") 
s1

## 5,6,11
s1 <- ggplot(data = d[d$variable == names(data)[11],], aes(x = factor(cl), y = value, fill = variable))
s1 <- s1 + geom_boxplot(fill = c('#FFCC00','#FF9900')) + geom_point(position = position_jitterdodge()) + stat_compare_means(method = "wilcox.test") 
s1



变量相关性
cor_ <- cor(data, method = 'pearson')
View(cor_)
write.table(symnum(cor_), file.choose(), col.names = T, sep = "\t")



cl <-read.table(file.choose(), header = T, sep = "\t")
cl$cl<-as.factor(cl$cl)

library(ggplot2)
ggplot(cl, aes(x=cl, y=年龄)) + geom_boxplot(outlier.colour="red", outlier.shape=4, outlier.size=2) + geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.2)

ggplot(cl, aes(x=cl, y=身高)) + geom_boxplot(outlier.colour="red", outlier.shape=4, outlier.size=2) + geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.2)

ggplot(cl, aes(x=cl, y=体重)) + geom_boxplot(outlier.colour="red", outlier.shape=4, outlier.size=2) + geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.2)

ggplot(cl, aes(x=cl, y=腰围)) + geom_boxplot(outlier.colour="red", outlier.shape=4, outlier.size=2) + geom_dotplot(binaxis='y', stackdir='center', stackratio=1.5, dotsize=0.2)


myboxplot <- function(x, data, col = NULL, xlab, pvalue="auto") {
    boxplot(x, data, axes = FALSE, col = col)
    axis(1, at = 1:2, labels =FALSE)
    text(1:2, y=par()$usr[3]-0.08*(par()$usr[4]-par()$usr[3]),
         srt=60, xpd=T, adj=1, labels = xlab)
    if (pvalue == "auto") {
        pvalue <- round(t.test(x, data=data)$p.value, 3)
    }

    if (!is.null(pvalue)) {
        plab <- paste("p =", pvalue)
        text(1.5, y = par()$usr[4]*1.05, xpd=T, label=plab, col=col)
    }
}
