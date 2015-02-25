###########################################################################################################
###########################################################################################################
#Script for the analyses of DisperMyz data
###########################################################################################################
###########################################################################################################


#Loading the datafile into R
#First you need to set the right working directory
getwd()
setwd("~/work/Rfichiers/Githuber/data")

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

#two clone-corrected datasets are build on the complete dataset. The first one is a 'conservative' clone 
#correction since we only keep on MLG per site
BRAcccons<-BRA[BRA$cc_cons==1,]
sum(table(BRAcccons$pop_ID)) #only 259 individuals left
#the second clone corrected dataset is less conservative, we only removed over-represented multicopies 
#MLG (see Mat & Meth for details)
BRAcc<-BRA[BRA$cc==1,]
sum(table(BRAcc$pop_ID)) #only 338 individuals left


#the analysis here are performed for the clone-corrected dataset, but the same analysis can 
#be performed with the other dataset by replacing BRAcc by the complete (BRA) or the 
#conservative clone-corrected dataset (BRAcccons) in the following line code, and then rerun 
#other lines of code
BRAt<-BRAcc #name of the input file

#converting data to a genind format
BRADE<-df2genind(BRAt[,14:27],ncode=3,ind.names=BRAt$sample_ID, 
                 pop=BRAt$pop_ID,missing=0,ploidy=1)
BRADE@other$xy<-BRAt[,4:5]
#determination of the number of clusters
clustBRADE<- find.clusters(BRADE,max.n.clust=35)
#with 40 PCs, we lost nearly no information
clustBRADE<- find.clusters(BRADE,n.pca=40,max.n.clust=35) #chose 3 clusters
#which individuals in which clusters per population
table(pop(BRADE),clustBRADE$grp)
#DAPC by itself, first we try to optimized the number of principal component (PCs) 
#to retain to perform the analysis
dapcBRADE<-dapc(BRADE,clustBRADE$grp,n.da=5,n.pca=100)
temp<-optim.a.score(dapcBRADE)
dapcBRADE<-dapc(BRADE,clustBRADE$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcBRADE) #based on this result, we finaly chose 7 PCs
dapcBRADE<-dapc(BRADE,clustBRADE$grp,n.da=7,n.pca=7)
#STRUCTURE-like graphic
compoplot(dapcBRADE,lab=NA)
scatter(dapcBRADE,xax=1, yax=2)

BRADEpop<-genind2genpop(BRADE,process.other=T,missing="0")

image(alt,col=brewer.pal(9,"Greys"))
stars(table(pop(BRADE),dapcBRADE$assign),draw.segment=TRUE,
      locations=BRADEpop@other$xy,
      #locations=cbind(jitter(BRADEpop@other$xy$longitude,200),
      #                jitter(BRADEpop@other$xy$latitude,200)),
      add=T,len=0.5)




