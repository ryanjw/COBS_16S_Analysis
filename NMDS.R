fracs<-read.table("~/Google Drive/COBS_16S_biodiversity/phylo_and_fracs_unrar.txt",header=TRUE,check.names=FALSE,sep=" ")
# fracs<-read.table("phylo_and_fracs_unrar.txt",sep=" ",header=T,check.names=F)
ws<-read.table("~/Google Drive/COBS_16S_biodiversity/ws_unrared_table_singletons_removed.txt",header=TRUE,check.names=FALSE,sep=" ")
meta<-read.csv("~/Google Drive/COBS_16S_biodiversity/KBase_MGRast_Metadata_9May2013_EMB.csv",header=T)
names(meta)[1]<-"Sample"
# ws<-read.table("ws_unrared_table_singletons_removed.txt",header=TRUE,check.names=FALSE,sep=" ")
library(reshape2)
library(vegan)
library(lme4)
library(lmerTest)
ws_only<-subset(ws, SoilFrac=="WS")

#min(subset(ws_only, value > 0)$value)

ws_only$value<-ws_only$value/1000

combined<-rbind(fracs,ws_only)

dataset<-dcast(combined, Sample+Date+Block+Crop+SoilFrac~variable, value.var="value",fun.aggregate=sum)
rared<-rrarefy(dataset[,-c(1:5)],10903)
mds<-metaMDS(decostand(rared,"total"),k=2,autotransform=F)
mds
MDS1<-data.frame(scores(mds))$NMDS1
MDS2<-data.frame(scores(mds))$NMDS2
Treatment<-dataset$SoilFrac
Crop<-dataset$Crop
NMDS<-data.frame(MDS1,MDS2,Treatment,Crop)

# function for making  an ellipse
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 1000) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }

  df_ell <- data.frame()
  for(g in levels(NMDS$Treatment)){
    df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$Treatment==g,],
                    veganCovEllipse(cov.wt(cbind(MDS1,MDS2),wt=rep(1/length(MDS1),length(MDS1)))$cov,center=c(mean(MDS1),mean(MDS2)))))
                    ,group=g))
  }
head(df_ell)

df_ell$group<-factor(df_ell$group, levels=c("Micro","SM","MM","LM","WS"))
levels(df_ell$group)<-c("Micro","Small","Medium","Large","Whole Soil")
NMDS$Treatment<-factor(NMDS$Treatment, levels=c("Micro","SM","MM","LM","WS"))
levels(NMDS$Treatment)<-c("Micro","Small","Medium","Large","Whole Soil")
NMDS$Crop<-factor(NMDS$Crop)
levels(NMDS$Crop)<-c("Corn","Prairie","Fertilized\nPrairie")

ggplot(data = NMDS, aes(MDS1, MDS2)) + geom_point(aes(color = Treatment,shape=Crop),size=3,alpha=0.75)+geom_path(data=df_ell, aes(x=MDS1, y=MDS2,colour=group), size=1, linetype=5)+scale_color_manual(name="Aggregate\nFraction",values=(c("black","#d6604d","#b2182b","#67001f","#878787")))+theme_bw(base_size=17)+theme(aspect.ratio=1)+annotate("text",x=.7,y=.25, label="Stress = 0.097")
