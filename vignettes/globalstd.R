## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----demodata-----------------------------------------------------------------
library(pmd)
data("spmeinvivo")
str(spmeinvivo)


## -----------------------------------------------------------------------------
knitr::include_graphics('https://yufree.github.io/presentation/figure/GlobalStd.png')


## ----rtg----------------------------------------------------------------------
pmd <- getpaired(spmeinvivo)
plotrtg(pmd)


## ----pmd----------------------------------------------------------------------
plotpaired(pmd)


## ----pmdindex-----------------------------------------------------------------
# show the unique PMD found by getpaired function
for(i in 1:length(unique(round(pmd$paired$diff,2)))){
        diff <- unique(round(pmd$paired$diff,2))[i]
        index <- round(pmd$paired$diff,2)== diff
        plotpaired(pmd,index)
}


## ----std----------------------------------------------------------------------
std <- getstd(pmd)


## ----stdplot------------------------------------------------------------------
plotstd(std)


## ----stdrtplot----------------------------------------------------------------
par(mfrow = c(2,3))
plotstdrt(std,rtcluster = 23,main = 'Retention time group 23')
plotstdrt(std,rtcluster = 9,main = 'Retention time group 9')
plotstdrt(std,rtcluster = 18,main = 'Retention time group 18')
plotstdrt(std,rtcluster = 67,main = 'Retention time group 67')
plotstdrt(std,rtcluster = 49,main = 'Retention time group 49')
plotstdrt(std,rtcluster = 6,main = 'Retention time group 6')


## ----globalcor----------------------------------------------------------------
std2 <- globalstd(pmd,corcutoff = 0.9)


## ----pca----------------------------------------------------------------------
library(enviGCMS)
par(mfrow = c(2,2),mar = c(4,4,2,1)+0.1)
plotpca(std$data,lv = as.numeric(as.factor(std$group$sample_group)),main = "all peaks")
plotpca(std$data[std$stdmassindex,],lv = as.numeric(as.factor(std$group$sample_group)),main = paste(sum(std$stdmassindex),"independent peaks"))
plotpca(std2$data[std2$stdmassindex,],lv = as.numeric(as.factor(std$group$sample_group)),main = paste(sum(std2$stdmassindex),"reduced independent peaks"))


## ----comp---------------------------------------------------------------------
stdcluster <- getpseudospectrum(std)
# extract the first pseudospectra for retention time cluster 37
idx <- stdcluster$pseudo$sid=='37@1'
plot(stdcluster$pseudo$mz[idx],stdcluster$pseudo$ins[idx],type = 'h',xlab = 'm/z',ylab = 'intensity',main = 'pseudo spectra for retention time cluster 37')
# considering the correlation coefficient cutoff
stdcluster2 <- getpseudospectrum(std, corcutoff = 0.9)


## ----corpeak------------------------------------------------------------------
corcluster <- getcorpseudospectrum(spmeinvivo)
# extract pseudospectra 37@1
peak <- corcluster$pseudo[corcluster$pseudo$sid == '37@1',]
plot(peak$ins~peak$mz,type = 'h',xlab = 'm/z',ylab = 'intensity',main = 'pseudo spectra for correlation cluster')


## ----corcomp------------------------------------------------------------------
par(mfrow = c(2,3),mar = c(4,4,2,1)+0.1)
plotpca(std$data[std$stdmassindex,],lv = as.numeric(as.factor(std$group$sample_group)),main = paste(sum(std$stdmassindex),"independent peaks"))
plotpca(std$data[stdcluster$stdmassindex2,],lv = as.numeric(as.factor(std$group$sample_group)),main = paste(sum(stdcluster$stdmassindex2),"independent base peaks"))
plotpca(std$data[stdcluster2$stdmassindex2,],lv = as.numeric(as.factor(std$group$sample_group)),main = paste(sum(stdcluster2$stdmassindex2),"independent reduced base peaks"))
plotpca(std$data[corcluster$stdmassindex,],lv = as.numeric(as.factor(std$group$sample_group)),main = paste(sum(corcluster$stdmassindex),"peaks without correlationship"))
plotpca(std$data[corcluster$stdmassindex2,],lv = as.numeric(as.factor(std$group$sample_group)),main = paste(sum(corcluster$stdmassindex2),"base peaks without correlationship"))
plotpca(std$data,lv = as.numeric(as.factor(std$group$sample_group)),main = paste(nrow(std$data),"all peaks"))

# Evaluate multiple algorithms against original data
getsim(std$data, std$data[std$stdmassindex,])
getsim(std$data, std$data[stdcluster$stdmassindex2,])
getsim(std$data, std$data[corcluster$stdmassindex2,])


## ----nse----------------------------------------------------------------------
# Calculate SCM/NSE for the full data and the reduced standard masses
getnse(std$data)
getnse(std$data[std$stdmassindex,])


## ----sda----------------------------------------------------------------------
sda <- getsda(std)


## -----------------------------------------------------------------------------
library(igraph)
cdf <- sda$sda
# get the PMDs and frequency
pmds <- as.numeric(names(sort(table(cdf$diff2),decreasing = T)))
freq <- sort(table(cdf$diff2),decreasing = T)
# filter the frequency larger than 10 for demo
pmds <- pmds[freq>10]
cdf <- sda$sda[sda$sda$diff2 %in% pmds,]
g <- igraph::graph_from_data_frame(cdf,directed = F)
l <- igraph::layout_with_fr(g)
for(i in 1:length(pmds)){
  g2 <- igraph::delete_edges(g,which(E(g)$diff2%in%pmds[1:i]))
  plot(g2,edge.width=1,vertex.label="",vertex.size=1,layout=l,main=paste('Top',length(pmds)-i,'high frequency PMDs'))
}


## ----stdsda-------------------------------------------------------------------
plotstdsda(sda)


## ----stdsdaidx----------------------------------------------------------------
par(mfrow = c(1,3),mar = c(4,4,2,1)+0.1)
plotstdsda(sda,sda$sda$diff2 == 2.02)
plotstdsda(sda,sda$sda$diff2 == 28.03)
plotstdsda(sda,sda$sda$diff2 == 58.04)


## ----all----------------------------------------------------------------------
sdaall <- getsda(spmeinvivo)
par(mfrow = c(1,3),mar = c(4,4,2,1)+0.1)
plotstdsda(sdaall,sdaall$sda$diff2 == 2.02)
plotstdsda(sdaall,sdaall$sda$diff2 == 28.03)
plotstdsda(sdaall,sdaall$sda$diff2 == 58.04)


## -----------------------------------------------------------------------------
sda2 <- getsda(std, corcutoff = 0.9)
plotstdsda(sda2)


## ----rda----------------------------------------------------------------------
sda <- getrda(spmeinvivo$mz)
# check high frequency pmd
colnames(sda)
# get certain pmd related m/z
idx <- sda[,'2.016']
# show the m/z
spmeinvivo$mz[idx]


## ----wrap---------------------------------------------------------------------
result <- globalstd(spmeinvivo, sda=FALSE)


## ----gettarget----------------------------------------------------------------
# you need retention time for independent peaks
index <- gettarget(std$rt[std$stdmassindex])
# output the ions for each injection
table(index)
# show the ions for the first injection
std$mz[index==1]
std$rt[index==1]

