## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup--------------------------------------------------------------------
library(pmd)
data("spmeinvivo")


## ----tarnet-------------------------------------------------------------------
library(igraph)
# check metabolites of C18H39NO
# Use common PMDs for biological reactions
chain <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99,58.04,13.98),mass = 286.3101,digits = 2,corcutoff = 0)
# show as network
net <- graph_from_data_frame(chain$sdac,directed = F)
pal <- grDevices::rainbow(5)
plot(net,vertex.label=round(as.numeric(V(net)$name),2),vertex.size =5,edge.width = 3,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],vertex.label.dist=1,vertex.color=ifelse(round(as.numeric(V(net)$name),4) %in% 286.3101,'red','black'), main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)
# Consider the correlation coefficient cutoff
chain <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99,58.04,13.98),mass = 286.3101,digits = 2,corcutoff = 0.6)
# show as network
net <- graph_from_data_frame(chain$sdac,directed = F)
pal <- grDevices::rainbow(5)
plot(net,vertex.label=round(as.numeric(V(net)$name),2),vertex.size =5,edge.width = 3,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],vertex.label.dist=1,vertex.color=ifelse(round(as.numeric(V(net)$name),4) %in% 286.3101,'red','black'), main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)


## ----net----------------------------------------------------------------------
std <- globalstd(spmeinvivo,sda = F)
sda <- getsda(std,freqcutoff = 12)
df <- sda$sda
net <- graph_from_data_frame(df,directed = F)
pal <- grDevices::rainbow(length(unique(E(net)$diff2)))
plot(net,vertex.label=NA,vertex.size = 5,edge.width = 3,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)


## ----nwa----------------------------------------------------------------------
# network community structure detection
ceb <- cluster_edge_betweenness(net,weights = abs(E(net)$cor), directed = F) 
plot(ceb, net,vertex.label=NA,vertex.size = 5,edge.width = 3,) 
# output membership
head(cbind(ceb$membership,ceb$names))


## ----cda----------------------------------------------------------------------
cbp <- enviGCMS::getfilter(std,rowindex = std$stdmassindex)
cda <- getcda(cbp)
df <- cda$cda
# filter based on retention time differences larger than 2 mins
df <- df[df$diffrt>120,]
netc <- graph_from_data_frame(df,directed = F)
plot(netc,vertex.label=NA,vertex.size = 5,edge.width = 3,main = 'Correlation network')


## ----chainseq-----------------------------------------------------------------
# Trace homologous series with 14.01 PMD (+CH2)
homo <- gethomolog(spmeinvivo, unit = 14.01)

# Trace complex sequences e.g. oxidation (15.99) followed by multiple repeated methylation (14.01+)
sequence <- getchainseq(spmeinvivo, "15.99; 14.01+")


## ----source-------------------------------------------------------------------
deg <- degree(net, mode = 'all')
median(deg)
endogenous <- names(deg)[deg>median(deg)]
exogenous <- names(deg)[deg<=median(deg)]


## -----------------------------------------------------------------------------
pmd <- getreact(spmeinvivo,pmd=15.99)
# show the ions with the same PMD
head(pmd$pmd)
# show the corresponding quantitative PMD data across samples, each row show the sum of intensity of paired masses qualified for stable mass pairs
head(pmd$pmddata)


## -----------------------------------------------------------------------------
spmeinvivo$rt <- NULL
pmd <- getreact(spmeinvivo,pmd=15.99)
# show the ions with the same PMD
head(pmd$pmd)
# show the corresponding quantitative PMD data across samples, each row show the sum of intensity of paired masses qualified for stable mass pairs
head(pmd$pmddata)


## -----------------------------------------------------------------------------
data("spmeinvivo")
pmd <- getreact(spmeinvivo,pmd=15.99,method = 'dynamic')
# show the ions with the same PMD
head(pmd$pmd)
# show the corresponding quantitative PMD data across samples, each row show the sum of intensity of paired masses qualified for stable mass pairs
head(pmd$pmddata)


## -----------------------------------------------------------------------------
data("spmeinvivo")
# remove redundant peaks
list <- globalstd(spmeinvivo,sda = T)
newlist <- enviGCMS::getfilter(list,rowindex = list$stdmassindex)
# get high frequency pmd
hfpmd <- unique(newlist$sda$diff2)
# generate quantitative results
pmd <- getreact(newlist,pmd=hfpmd)
# output the kegg pmd in the data
table(pmd$pmd$diff2)
# output quantitative result for each PMD
head(pmd$pmddata)
# output quantitative result for unique PMD
upmd <- aggregate(pmd$pmddata, by=list(pmd$pmd$diff2),sum)
# column for samples and row for unique PMD
head(upmd)


## -----------------------------------------------------------------------------
# output all existing PMD in KEGG
keggpmd <- unique(round(keggrall$pmd,2))
data("spmeinvivo")
# remove redundant peaks
list <- globalstd(spmeinvivo)
newlist <- enviGCMS::getfilter(list,rowindex = list$stdmassindex)
# generate quantitative results
pmd <- getreact(newlist,pmd=keggpmd)
# output the kegg pmd in the data
table(pmd$pmd$diff2)
# output quantitative result for each PMD
head(pmd$pmddata)
# output quantitative result for unique PMD
upmd <- aggregate(pmd$pmddata, by=list(pmd$pmd$diff2),sum)
# column for samples and row for unique PMD
head(upmd)


## -----------------------------------------------------------------------------
data(spmeinvivo)
# get the m/z
mz <- spmeinvivo$mz
# get the m/z intensity for all m/z, the row order is the same with mz
insms <- spmeinvivo$data
# check high frequency pmd
sda <- getrda(mz)
colnames(sda)
# save them as numeric vector
hfpmd <- as.numeric(colnames(sda))


## -----------------------------------------------------------------------------
# get details for certain pmd
pmddf <- getpmddf(mz,pmd=18.011,digits = 3)
# add intensity for all the paired ions
mz1ins <- insms[match(pmddf$ms1,mz),]
mz2ins <- insms[match(pmddf$ms2,mz),]
# get the pmd pair intensity
pmdins <- mz1ins+mz2ins
# get the pmd total intensity across samples
pmdinsall <- apply(pmdins,2,sum)
# show the PMD intensity
pmdinsall


## -----------------------------------------------------------------------------
# get the ratio of larger m/z over smaller m/z
ratio <- mz2ins/mz1ins
# filter PMD based on RSD% across samples
# cutoff 30%
cutoff <- 0.3
# get index for static PMD
rsdidx <- apply(ratio,1,function(x) sd(x)/mean(x)<cutoff)
# get static PMD
pmddfstatic <- pmddf[rsdidx,]
# get static intensity
pmdinsstatic <- pmdins[rsdidx,]
# normalize the ions pair intensity to avoid influences from large response factors
pmdinsstaticscale <- t(scale(t(pmdinsstatic)))
# get the pmd static intensity across samples
pmdinsstaticall <- apply(pmdinsstaticscale,2,sum)
# show the PMD static intensity for each sample
pmdinsstaticall

# get index for dynamic PMD
rsdidx <- apply(ratio,1,function(x) sd(x)/mean(x)>=cutoff)
# get dynamic PMD
pmddfdynamic <- pmddf[rsdidx,]
# get dynamic intensity for ms1 and ms2
pmdinsdynamicms1 <- apply(mz1ins[rsdidx,],1,function(x) sd(x)/mean(x))
pmdinsdynamicms2 <- apply(mz2ins[rsdidx,],1,function(x) sd(x)/mean(x))
# find the stable ms and use ratio as intensity
idx <- pmdinsdynamicms1>pmdinsdynamicms2
pmdinsdynamic <- ratio[rsdidx,]
pmdinsdynamic[idx,] <- 1/ratio[rsdidx,][idx,]
# get the pmd dynamic intensity across samples
pmdinsdynamicall <- apply(pmdinsdynamic,2,sum)
# show the PMD dynamic intensity for each sample
pmdinsdynamicall


## -----------------------------------------------------------------------------
# get details for certain pmd
pmddf <- getpmddf(mz,pmd=hfpmd,digits = 3)
# viz by igraph package
library(igraph)
net <- graph_from_data_frame(pmddf,directed = F)
pal <- grDevices::rainbow(length(unique(E(net)$diff2)))
plot(net,vertex.label=NA,vertex.size = 5,edge.width = 3,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)


## -----------------------------------------------------------------------------
data(spmeinvivo)
spmeinvivo$rt <- NULL
chain <- getchain(spmeinvivo,diff = c(2.02,14.02,15.99,58.04,13.98),mass = 286.3101,digits = 2,corcutoff = 0)
# show as network
net <- graph_from_data_frame(chain$sdac,directed = F)
pal <- grDevices::rainbow(5)
plot(net,vertex.label=round(as.numeric(V(net)$name),2),vertex.size =5,edge.width = 3,edge.color = pal[as.numeric(as.factor(E(net)$diff2))],vertex.label.dist=1,vertex.color=ifelse(round(as.numeric(V(net)$name),4) %in% 286.3101,'red','black'), main = 'PMD network')
legend("topright",bty = "n",
       legend=unique(E(net)$diff2),
       fill=unique(pal[as.numeric(as.factor(E(net)$diff2))]), border=NA,horiz = F)


## -----------------------------------------------------------------------------
# all reaction
data("omics")
head(omics)
# kegg reaction
data("keggrall")
head(keggrall)
# literature reaction for mass spectrometry
data("sda")
head(sda)


## -----------------------------------------------------------------------------
data("hmdb")
head(hmdb)


## -----------------------------------------------------------------------------
plotcn('C6H12O6','Glucose',c(2.016,14.016,15.995))

