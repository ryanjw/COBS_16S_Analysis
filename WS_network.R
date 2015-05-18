library(reshape2)
library(vegan)
library(plyr)
ws<-read.table("ws_unrared_table_singletons_removed.txt",header=T,sep=" ")
meta<-read.csv("KBase_MGRast_Metadata_9May2013_EMB.csv")
names(meta)[1]<-"Sample"
ws$value<-ws$value/1000
ws<-subset(ws, SoilFrac=="WS")
frac_df<-dcast(ws, Sample+Date+Block+Crop+SoilFrac~genus, value.var="value",fun.aggregate=sum)

results<-data.frame()

options(warnings=-1)
for(a in 1:100){
	counter<-0
	frac_df_rar<-cbind(frac_df[,1:5],rrarefy(frac_df[,-c(1:5)],10903))
	frac_merge<-merge(frac_df_rar,meta[,c(1,6,8)],by="Sample")
	dataset<-frac_merge
	

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
				new.row<-data.frame(a,"WS",names(dataset)[b],names(dataset)[c],rho,p.value,species1.ab,species2.ab)
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
