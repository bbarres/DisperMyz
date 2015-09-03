###############################################################################
###############################################################################
#R code for FAM report figures
###############################################################################
###############################################################################

#loading the packages necessary for the analysis
library(adegenet)
library(gdata)
library(RColorBrewer)

#Setting the right working directory
setwd("~/work/Rfichiers/Githuber/DisperMyz_data")


###############################################################################
#loading and preparing the dataset
###############################################################################

#first, we load the genetic dataset. It consists in the unique MLG for the 
#entire dataset
MyzFAM<-read.table("FAM.dat",header=T,sep="\t")
#here is the structure of the datafile, for explanation of each columns, see 
#ReadMe.txt file in DRYAD repository
head(MyzFAM)
#a summary of the different variables
summary(MyzFAM)
colnames(MyzFAM)
#total number of individuals
dim(MyzFAM)[1] #682 individuals

JDD<-MyzFAM #name of the input file
JDD<-drop.levels(JDD)
#let's define a set of color for keeping some consistency in the plots
coloor<-c("orange","green","blue","yellow","hotpink")


###############################################################################
#DAPC on microsatellites only
###############################################################################

#converting data to a genind format, first we use only the microsatellite data
JDDmicro<-df2genind(JDD[,c("MP_27","MP_39","MP_44","MP_5","MP_7","MP_23",
                           "MP_45","MP_28","MP_9","MP_13","MP_2","MP_38",
                           "MP_4","MP_46")],
                    ncode=6,ind.names=JDD$individual, 
                    pop=JDD$host,ploidy=2)
#include the coordinates of the samples
JDDmicro@other$xy<-JDD[,c("longitude","latitude")]
#we can also include the resistance genotypes as supplementary information
JDDmicro@other$KDR<-JDD[,"KDR"]
JDDmicro@other$sKDR<-JDD[,"sKDR"]
JDDmicro@other$MACE<-JDD[,"MACE"]
JDDmicro@other$R81T<-JDD[,"R81T"]

#now we analyse the adegenet format dataset with dapc
JDDade<-JDDmicro
#determination of the number of clusters
clustJDDade<-find.clusters(JDDade,max.n.clust=30)
#with 50 PCs, we lost nearly no information and after K=5, the decrease of 
#the BIC value is smaller, so we chose the maximum number of clusters to be 5 
#which individuals in which clusters per population
table(pop(JDDade),clustJDDade$grp)
#We try to optimize the number of principal component (PCs) to retain to 
#perform the analysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=5,n.pca=30)
temp<-optim.a.score(dapcJDDade)
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=4,n.pca=15)
temp<-optim.a.score(dapcJDDade)
#the optimal number of PCs fell between 5 and 9 (depending on the run), so we 
#chose the smallest number of PCs (5) in order to avoid overfitting of the 
#model. Then we do the actual DAPC anlysis
dapcJDDade<-dapc(JDDade,clustJDDade$grp,n.da=3,n.pca=10)
#STRUCTURE-like graphic
compoplot(dapcJDDade,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
#scatter plot
scatter(dapcJDDade,xax=1, yax=2,col=coloor)
scatter(dapcJDDade,xax=1, yax=3,col=coloor)
scatter(dapcJDDade,xax=2, yax=3,col=coloor)

#Run the 'find.clusters' and DAPC analysis for K=2 to 5
clustJDDade2<-find.clusters(JDDade,n.pca=50,n.clust=2)
dapcJDDade2<-dapc(JDDade,clustJDDade2$grp,n.da=1,n.pca=2)
compoplot(dapcJDDade2,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
clustJDDade3<-find.clusters(JDDade,n.pca=50,n.clust=3)
dapcJDDade3<-dapc(JDDade,clustJDDade3$grp,n.da=2,n.pca=3)
compoplot(dapcJDDade3,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
clustJDDade4<-find.clusters(JDDade,n.pca=50,n.clust=4)
dapcJDDade4<-dapc(JDDade,clustJDDade4$grp,n.da=3,n.pca=4)
compoplot(dapcJDDade4,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)
clustJDDade5<-find.clusters(JDDade,n.pca=50,n.clust=5)
dapcJDDade5<-dapc(JDDade,clustJDDade5$grp,n.da=3,n.pca=5)
compoplot(dapcJDDade5,lab=pop(JDDade),legend=FALSE,
          cex.names=0.3,cex.lab=0.5,cex.axis=0.5,col=coloor)

#a more beautifull scatter plot, the colors are matching colors used in 
#the structure-like plot
scatter(dapcJDDade5,xax=1,yax=2,cstar=1,cell=0,clab=0,main="Axis 1 & 2",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)
scatter(dapcJDDade5,xax=2,yax=3,cstar=1,cell=0,clab=0,main="Axis 2 & 3",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)

#the same plot with resistance genotypes added
scatter(dapcJDDade5,xax=1,yax=2,cstar=1,cell=0,clab=0,main="Axis 1 & 2",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)
#adding KDR resistotype
points(dapcJDDade5$ind.coord[,1],dapcJDDade5$ind.coord[,2],pch=21,xpd=NA,
       col="black",cex=1.5,bg=as.numeric(as.factor(JDDmicro@other$KDR)))
legend("topright",levels(as.factor(JDDmicro@other$KDR)),col="black",pch=21,
       pt.bg=levels(as.factor(as.numeric(as.factor(JDDmicro@other$KDR)))),
       xpd=NA)
title("KDR resistotypes")

#adding sKDR resistotype
scatter(dapcJDDade5,xax=1,yax=2,cstar=1,cell=0,clab=0,main="Axis 1 & 2",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade5$ind.coord[,1],dapcJDDade5$ind.coord[,2],pch=21,xpd=NA,
       col="black",cex=1.5,bg=as.numeric(as.factor(JDDmicro@other$sKDR)))
legend("topright",levels(as.factor(JDDmicro@other$sKDR)),col="black",pch=21,
       pt.bg=levels(as.factor(as.numeric(as.factor(JDDmicro@other$sKDR)))),
       xpd=NA)
title("sKDR resistotypes")

#adding MACE resistotype
scatter(dapcJDDade5,xax=1,yax=2,cstar=1,cell=0,clab=0,main="Axis 1 & 2",
        solid=0.3,col=rainbow(5)[c(2,3,5,1,4)],pch=19,cex=3,scree.da=FALSE)
points(dapcJDDade5$ind.coord[,1],dapcJDDade5$ind.coord[,2],pch=21,xpd=NA,
       col="black",cex=1.5,bg=as.numeric(as.factor(JDDmicro@other$MACE)))
legend("topright",levels(as.factor(JDDmicro@other$MACE)),col="black",pch=21,
       pt.bg=levels(as.factor(as.numeric(as.factor(JDDmicro@other$MACE)))),
       xpd=NA)
title("MACE resistotypes")

#in case we need the q-matrix of the individuals for other purposes
write.table(dapcJDDade2$posterior,file="AgrAphDAPCK2.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)
write.table(dapcJDDade3$posterior,file="AgrAphDAPCK3.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)
write.table(dapcJDDade4$posterior,file="AgrAphDAPCK4.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)
write.table(dapcJDDade5$posterior,file="AgrAphDAPCK5.txt",sep="\t",
            quote=FALSE,row.names=TRUE,col.names=FALSE)


###############################################################################
#Defining a function to make structure-like plot
###############################################################################

structplot<-function(qmat,coolcol,effP,nameP,leg_y="",cexpop=1,cexy=2,
                     mef=c(1,1,1,1,1),colbord="NA",angl=0,distxax=0.005)
  #'qmat': the q-matrix like matrix
  #'coolcol': a vector of colors, for the different genetic clusters
  #'effP': a vector giving the number of individuals in each population
  #'nameP': a list giving the names of the different populations
  #'leg_y': a characters string used for the Y-axis legend
  #'cexpop': the cex factor for the population name on the x-axis
  #'cexy': the cex factor for the Y-legend
  #'mef': a vector of length 5 to pimp the graph. Each 1 value add a feature
  #the first is for the external rectangle, the second is to deliminate the 
  #different populations, the third is for adding an x-axis with tick, the 
  #fourth is for adding the name of the different populations and the fifth 
#is for adding a Y legend
#'colbord': the color of the line between individuals
#'angl': the angle of the tag of the x-axis
#'distxax': control the distance of the tag to the x-axis

{
  barplot(qmat,col=coolcol,beside=FALSE,border=colbord,
          space=0,ylim=c(-0.03,1.03),axisnames=FALSE,axes=FALSE)
  
  if(mef[1]==1) {
    #drawing an external rectangle
    rect(0-1/dim(qmat)[2],
         0-1/500,
         dim(qmat)[2]+1/dim(qmat)[2],
         1+1/500,
         lwd=3)
  }
  
  if(mef[2]==1) {
    #deliminated the different populations
    rect(c(0,cumsum(effP))[1:length(effP)],
         rep(0,length(effP)),
         cumsum(effP),
         rep(1,length(effP)),
         lwd=2)
  }
  
  if(mef[3]==1) {
    #add an x-axis
    axis(1,at=c(0,cumsum(effP))[1:length(effP)]+
           (cumsum(effP)-c(0,cumsum(effP))[1:length(effP)])/2,
         labels=FALSE,pos=0,lwd.ticks=2)
  }
  
  if(mef[4]==1) {
    #add the name of the different populations
    text(c(0,cumsum(effP))[1:length(effP)]+
           (cumsum(effP)-c(0,cumsum(effP))[1:length(effP)])/2,
         rep(par("usr")[3]-distxax,length(effP)),labels=nameP,srt=angl,
         xpd=NA,pos=1,cex=cexpop)
  }
  
  if(mef[5]==1) {
    #add some legend on the Y-axis
    mtext(leg_y,side=2,las=1,cex=cexy,adj=0.5,line=1)
  }
  
}

#some examples of the use of the function
#first you need to gather the number of individuals in each populations
effpop<-c(380,71,28,17,168,18)
#the names of the different populations might be useful too
poptiquet<-c("Pêcher","Colza","Tabac","Autres\nSecondaires","Piège\naérien",
             "Hôtes\nmultiples")
#be careful to use the same dataset that has been used for the DAPC 
#computation
structplot(t(dapcJDDade5$posterior),rainbow(5),effpop,poptiquet)
structplot(t(dapcJDDade5$posterior),rainbow(5),effpop,poptiquet,
           colbord="grey70",leg_y="K=5",angl=30,distxax=0.01)
structplot(t(dapcJDDade5$posterior),coloor,effpop,poptiquet,
           leg_y="K=5",mef=c(1,0,1,1,1),cexpop=0.5,cexy=5)
structplot(t(dapcJDDade2$posterior),coloor,effpop,poptiquet,
           colbord=0,leg_y="K=2",mef=c(0,1,0,1,1))

#Now, we can easily plot several structure-like plot in the same figure
op<-par(mfrow=c(4,1),mar=c(0,4,0,0),oma=c(3,0,0,0))
structplot(t(dapcJDDade5$posterior)[c(4,1,2,5,3),],rainbow(5),effpop,poptiquet,
           leg_y="K=5",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(t(dapcJDDade4$posterior)[c(1,4,3,2),],rainbow(5),effpop,poptiquet,
           leg_y="K=4",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(t(dapcJDDade3$posterior)[c(1,3,2),],rainbow(5),effpop,poptiquet,
           leg_y="K=3",cexy=1.2,mef=c(0,1,0,0,1),colbord="grey70")
structplot(t(dapcJDDade2$posterior),rainbow(5),effpop,poptiquet,
           leg_y="K=2",cexy=1.2,mef=c(0,1,1,1,1),colbord="grey70",
           distxax=0.08)
par(op)

#export to pdf 15 X 4 inches

#alternatively, you can import a q-matrix file and use the function in the 
#same manner. Be careful howerer to respect the order of the individuals and 
#the order of their respective populations
strK2<-t(read.table("FAMoutK2.str",header=FALSE,sep="\t")[,c(-1)])
strK3<-t(read.table("FAMoutK3.str",header=FALSE,sep="\t")[,c(-1)])
strK4<-t(read.table("FAMoutK4.str",header=FALSE,sep="\t")[,c(-1)])

display.brewer.all()
coloor <- brewer.pal(4,"Greys")

coloor <- c("firebrick","forestgreen","dodgerblue3","khaki2")

op<-par(mfrow=c(3,1),mar=c(0,3,0,0),oma=c(8,0,0,0))
structplot(strK4,coloor,effpop,poptiquet,
           leg_y="K=4",cexy=1.2,mef=c(0,1,0,0,1),colbord=NA)
structplot(strK3,coloor,effpop,poptiquet,
           leg_y="K=3",cexy=1.2,mef=c(0,1,0,0,1),colbord=NA)
structplot(strK2,coloor,effpop,poptiquet,
           leg_y="K=2",cexy=1.2,mef=c(0,1,1,1,1),colbord=NA,
           distxax=0.15,angl=45,cexpop=2)
par(op)


###############################################################################
#Identifying the best K for STRUCTURE run
###############################################################################

#Analyzes were performed using STRUCTURE2.3.4 software, with a model allowing 
#admixture and correlation of allele frequencies. Each run consisted of a 
#burn-in period of 10.000 iterations followed by 50.000 simulations. One 
#hundred repetitions of each run were performed for K ranging from 1 to 10. 
#Before importing the file, replace white space in the column header names 
#with underscore, replace "?1" by "alpha", and remove double white spaces or 
#it will provoc importation problem or failure

resstr_cccons<-read.table(file="FAM_deltaK.txt", header=T,sep="\t",
                          blank.lines.skip=T)

#a function which compute delta K values, nb_K is the number of different K 
#considered, and nb_rep is the number of repetition of each K
chooseK<-function(str_out,nb_K,nb_rep) {
  datatable<-data.frame("K"=c(rep(1:nb_K,each=nb_rep)),"Ln(Pd)"=str_out[,4])
  Lprim<-c(rep("NA",nb_rep))
  for (i in ((nb_rep+1):(nb_K*nb_rep))) {
    Lprim<-c(Lprim,str_out[i,4]-str_out[i-nb_rep,4])
  }
  datatable<-data.frame(datatable,as.numeric(Lprim))
  Lsecond<-c(rep("NA",nb_rep))
  for (i in (((2*nb_rep)+1):(nb_K*nb_rep))) {
    Lsecond<-c(Lsecond,abs(datatable[i,3]-datatable[i-nb_rep,3]))
  }
  Lsecond<-c(Lsecond,rep("NA",nb_rep))
  datatable<-data.frame(datatable,as.numeric(Lsecond))
  reztable<-data.frame("K"=c(1:nb_K))
  meanL<-c()
  sdL<-c()
  for (i in (1:nb_K)) {
    meanL<-c(meanL,mean(datatable[datatable$K==i,2]))
    sdL<-c(sdL,sd(datatable[datatable$K==i,2]))
  }
  reztable<-data.frame(reztable,meanL,sdL)
  meanLprime<-c()
  sdLprime<-c()
  for (i in (1:nb_K)) {
    meanLprime<-c(meanLprime,mean(as.numeric(datatable[datatable$K==i,3])))
    sdLprime<-c(sdLprime,sd(datatable[datatable$K==i,3]))
  }
  reztable<-data.frame(reztable,meanLprime,sdLprime)
  meanLsecond<-c()
  sdLsecond<-c()
  for (i in (1:nb_K)) {
    meanLsecond<-c(meanLsecond,mean(as.numeric(datatable[datatable$K==i,4])))
    sdLsecond<-c(sdLsecond,sd(datatable[datatable$K==i,4]))
  }
  reztable<-data.frame(reztable,meanLsecond,sdLsecond)
  deltaK<-c()
  for (i in (1:nb_K)) {
    deltaK<-c(deltaK,reztable[reztable$K==i,6]/reztable[reztable$K==i,3])
  }
  reztable<-data.frame(reztable,deltaK)
  return(reztable)
}

deltastr_cccons<-chooseK(resstr_cccons,10,10)

#a function to plot variation of Delta K and Ln(P(X|K)) with K. 
plotdeltaK<-function(datadeltaK,nb_K,titre){
  #'datadeltak': the output file of 'chooseK' function
  #'nb_K': the number of different K considered
  #'titre': the title of the plot you want to be displayed
  op<-par(pty="s")
  plot(datadeltaK[1:(nb_K-2),8],type="b",pch=24,cex=2.5,lwd=4,lty=1,
       col="transparent",bg="white",bty="n",ann=F)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),8],type="b",pch=24,bty="n",xaxt="n",yaxt="n",
       ann=F,cex=2.5,lwd=4,lty=1)
  axis(side=1,at=seq(1,13,1),lwd=3,font.axis=2)
  axis(side=2,lwd=3,font.axis=2)
  title(ylab="Delta K",font.lab=2,cex.lab=1.5)
  par(new=TRUE)
  plot(datadeltaK[1:(nb_K-2),2],type="b",pch=22,cex=2.5,lwd=4,lty=2,
       col="grey50",bg="white",bty="n",xaxt="n",yaxt="n",ann=F)
  axis(side=4,lwd=3,font.axis=2,col="grey50")
  mtext("Ln(P(X|K))", side=4, line=4,font=2,cex=1,col="grey50")
  title(main=titre,xlab="K",font.lab=2,cex.lab=1.5,cex.main=2)
  par(op)
}


plotdeltaK(deltastr_cccons,10,
           "Détermination du meilleur K par méthode Evanno (n=682)")

#you can obtain the same figure as in the manuscript by exporting the plot to 
#png format, with a width of 2400 X 1100 pixels



###############################################################################
#END
###############################################################################