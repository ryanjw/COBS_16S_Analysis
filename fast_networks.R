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

frac_df<-dcast(combined, Sample+Date+Block+Crop+SoilFrac~genus, value.var="value",fun.aggregate=sum)
write.table(frac_df, "~/Google Drive/COBS_16S_biodiversity/WS_and_Fracs_and_meta_at_genus.txt",sep=" ",row.names=F)
head(meta)
results<-data.frame()

options(warnings=-1)
for(a in 1:100){
	counter<-0
	frac_df_rar<-cbind(frac_df[,1:5],rrarefy(frac_df[,-c(1:5)],10903))
	frac_merge<-merge(frac_df_rar,meta[,c(1,6,8)],by="Sample")
	dataset<-frac_merge
	dataset<-subset(dataset, SoilFrac=="LM")

for(b in 7:(dim(dataset)[2]-1)){
		#every species will be compared to every other species, so there has to be another loop that iterates down the rest of the columns
	for(c in (b+1):(dim(dataset)[2])){
			
			#summing the abundances of species of the columns that will be compared
			species1.ab<-sum(dataset[,b],na.rm=T)
			species2.ab<-sum(dataset[,c],na.rm=T)

			#if the column is all 0's no co-occurrence will be performed
			if(species1.ab/12 >=1 & species2.ab/12 >=1){
				test<-cor.test(dataset[,b],dataset[,c],method="spearman",na.action=na.rm)
				rho<-test$estimate
				p.value<-test$p.value
				new.row<-data.frame(a,"LM",names(dataset)[b],names(dataset)[c],rho,p.value,species1.ab,species2.ab)
				results<-rbind(results,new.row)
			}
				
		}
		counter<-counter+1
		if(counter > 10){
		print(paste(b/dim(dataset)[2]*100,"% Done of Trial ",a,sep=""))		
		counter<-0
		}
}	
}


# trying faster way
library(Hmisc)

options(warnings=-1)
co_occur_results<-data.frame()
fracs<-as.vector(unique(frac_df$SoilFrac))
options(warnings=-1)
for(a in 1:100){
	counter<-0
	frac_df_rar<-cbind(frac_df[,1:5],rrarefy(frac_df[,-c(1:5)],10903))
	frac_merge<-merge(frac_df_rar,meta[,c(1,6,8)],by="Sample")
	dataset<-frac_merge
	
	for(f in 1:length(fracs)){
		dataset_sub<-subset(dataset, SoilFrac==fracs[f])
		thing<-rcorr(as.matrix(dataset_sub[,-c(1:6)]))
		dfr<-data.frame(thing$r,check.names=F)
		dfr<-cbind(names(dfr),dfr)
		names(dfr)[1]<-"n1"
	
		dfP<-data.frame(thing$P,check.names=F)
		dfP<-cbind(names(dfP),dfP)
		names(dfP)[1]<-"n1"
		
		dfr_melt<-melt(dfr)
		dfP_melt<-melt(dfP)
		names(dfr_melt)[2:3]<-c("n2","r")
		names(dfP_melt)[2:3]<-c("n2","P")
		df_melt_total<-merge(dfr_melt,dfP_melt, by=c("n1","n2"))
		df_melt_total<-na.omit(df_melt_total)
		df_melt_total$qval<-fdrtool(df_melt_total$P,statistic="pvalue",plot=F,verbose=F)$qval
		df_melt_total$SoilFrac<-fracs[f]
		co_occur_results<-rbind(co_occur_results,df_melt_total)
	}
	counter<-counter+1
	if(counter > 10){
		print(a/100)
		counter<-0
	}
}
	