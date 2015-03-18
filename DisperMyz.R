###########################################################################################################
###########################################################################################################
#Script for the analyses of DisperMyz data
###########################################################################################################
###########################################################################################################


#loading the packages necessary for the analysis
library(adegenet)
library(gdata)

#Loading the datafile into R, first you need to set the right working directory
getwd()
setwd("~/work/Rfichiers/Githuber/DisperMyz_data")

#first of all, we load the genetic dataset
DIMY<-read.table("dispermyz.dat",header=T,sep="\t")
#here is the structure of the datafile, for explanation of each columns, see ReadMe.txt file in DRYAD 
#repository
head(DIMY)
#a summary of the different variables
summary(DIMY)
colnames(DIMY)
#number of individuals in each sampled populations
table(DIMY$patch)
#total number of individuals
sum(table(DIMY$patch)) #850 individuals

#keeping only one MLG copy for the entire population
DIMYcccons<-DIMY[DIMY$dup=="o",]
sum(table(DIMYcccons$patch)) #only 404 individuals left
#keeping only one MLG copy for each patch
DIMYcc<-DIMY[DIMY$dup!="dp",]
sum(table(DIMYcc$patch)) #only 473 individuals left

#here you can choose one the three previous datafile (DIMY, DIMYcccons, DIMYcc)
JDD<-DIMYcc #name of the input file
JDD<-drop.levels(JDD)

#converting data to a genind format
JDDade<-df2genind(JDD[,c("MP27","MP39","MP44","MP5","MP7","MP23","MP45","MP28","MP9","MP13",
                         "MP2","MP38","MP4","MP46")],ncode=6,ind.names=JDD$ID_simple, 
                 pop=JDD$patch,missing=NA,ploidy=2)
#include the coordinates of the samples
JDDade@other$xy<-JDD[,c("longitude","latitude")]
#include the host from which individuals were sampled from
JDDade@other$host<-JDD[,c("host")]
#determination of the number of clusters
clustJDDade<- find.clusters(JDDade,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustJDDade<- find.clusters(JDDade,n.pca=40,max.n.clust=35) #chose 4 clusters
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#DAPC by itself, first we try to optimized the number of principal component (PCs) 
#to retain to perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=100)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=7)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=NA)
scatter(dapcJDDade,xax=1, yax=2)

coloor<-c("red","green","blue","orange")
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=TRUE)

#scatter plot with the different K groups and then plotting the sampling host
scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.3,pch=19,cex=3,scree.da=FALSE)
#oilseed_rape
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==1,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==1,2],col="black",
       ,cex=2,bg="black",pch=21)
#peach
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==3,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==3,2],col="black",
       ,cex=2,bg="black",pch=21)
#tobacco
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==4,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==4,2],col="black",
       ,cex=2,bg="black",pch=21)
#other_crop
points(dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==2,1],
       dapcJDDade$ind.coord[as.numeric(as.factor(JDDade@other$host))==2,2],col="black",
       ,cex=2,bg="black",pch=21)


scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.0,pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade$ind.coord[,1],dapcJDDade$ind.coord[,2],col=coloor[dapcJDDade$assign],
       pch=(as.numeric(as.factor(JDDade@other$host))+20),cex=2)

scatter(dapcJDDade,xax=1,yax=2,cstar=1,cell=0,clab=0,col=coloor,
        solid=0.0,pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade$ind.coord[,1],dapcJDDade$ind.coord[,2],col=coloor[dapcJDDade$assign],
       pch=21,bg=coloor[(as.numeric(as.factor(JDDade@other$host)))])


plot(JDDade@other$xy,cex=3,col=dapcJDDade$assign,pch=as.numeric(as.factor(JDDade@other$host)))



# BRADEpop<-genind2genpop(BRADE,process.other=T,missing="0")
# 
# image(alt,col=brewer.pal(9,"Greys"))
# stars(table(pop(JDDade),dapcJDDade$assign),draw.segment=TRUE,
#       locations=JDDade@other$xy,
#       #locations=cbind(jitter(BRADEpop@other$xy$longitude,200),
#       #                jitter(BRADEpop@other$xy$latitude,200)),
#       add=T,len=0.5)




