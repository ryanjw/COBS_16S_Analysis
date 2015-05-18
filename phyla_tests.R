fracs<-read.table("~/Google Drive/COBS_16S_biodiversity/phylo_and_fracs_unrar.txt",header=TRUE,check.names=FALSE,sep=" ")
# fracs<-read.table("phylo_and_fracs_unrar.txt",sep=" ",header=T,check.names=F)
ws<-read.table("~/Google Drive/COBS_16S_biodiversity/ws_unrared_table_singletons_removed.txt",header=TRUE,check.names=FALSE,sep=" ")
# ws<-read.table("ws_unrared_table_singletons_removed.txt",header=TRUE,check.names=FALSE,sep=" ")
library(reshape2)
library(vegan)
library(lme4)
library(lmerTest)
ws_only<-subset(ws, SoilFrac=="WS")

#min(subset(ws_only, value > 0)$value)

ws_only$value<-ws_only$value/1000

combined<-rbind(fracs,ws_only)


dataset<-dcast(combined, Sample+Date+Block+Crop+SoilFrac~phylum, value.var="value",fun.aggregate=sum)
dataset_prot<-dcast(combined, Sample+Date+Block+Crop+SoilFrac~class,value.var="value",fun.aggregate=sum)

final<-data.frame()
for(t in 1:100){
dataset_rar<-cbind(dataset[,1:5],rrarefy(dataset[,-c(1:5)],10903))
dataset_prot_rar<-cbind(dataset_prot[,1:5],rrarefy(dataset_prot[,-c(1:5)],10903))
dataset_prot_rar<-dataset_prot_rar[,c(1:5,23,30,42,50)]
dataset_comb<-cbind(dataset_rar,dataset_prot_rar[,-c(1:5)])

dataset_comb_melt<-melt(dataset_comb, id=c("Sample","Date","Block","Crop","SoilFrac"))
dataset_comb_melt<-subset(dataset_comb_melt,  variable=="  p__Acidobacteria"| variable=="  p__Verrucomicrobia"| variable=="  p__Actinobacteria"| variable=="  p__Bacteroidetes"| variable=="  p__Crenarchaeota"| variable=="  p__Chloroflexi"| variable=="  p__Gemmatimonadetes"| variable=="  p__Planctomycetes"| variable=="  p__Firmicutes"| variable=="  p__Nitrospirae"| variable=="  c__Alphaproteobacteria"| variable=="  c__Betaproteobacteria"| variable=="  c__Deltaproteobacteria"| variable=="  c__Gammaproteobacteria")

phys<-as.vector(unique(dataset_comb_melt$variable))

for(i in 1:length(phys)){
	#i<-1
	test<-lmer((value)~Date*Crop*SoilFrac+(1|Date/Crop/SoilFrac)+(1|Block),data=subset(dataset_comb_melt, variable==phys[i]),REML=FALSE)
	table<-data.frame(anova(test))
	table$taxon<-phys[i]
	final<-rbind(final,table)
}
print(t/100)
}


# now producing the estimates for crop*soilfrac and soilfrac
library(plyr)
final_int<-data.frame()
for(t in 1:100){
dataset_rar<-cbind(dataset[,1:5],rrarefy(dataset[,-c(1:5)],10903))
dataset_prot_rar<-cbind(dataset_prot[,1:5],rrarefy(dataset_prot[,-c(1:5)],10903))
dataset_prot_rar<-dataset_prot_rar[,c(1:5,23,30,42,50)]
dataset_comb<-cbind(dataset_rar,dataset_prot_rar[,-c(1:5)])

dataset_comb_melt<-melt(dataset_comb, id=c("Sample","Date","Block","Crop","SoilFrac"))
dataset_comb_melt<-subset(dataset_comb_melt,  variable=="  p__Acidobacteria"| variable=="  p__Verrucomicrobia"| variable=="  p__Actinobacteria"| variable=="  p__Bacteroidetes"| variable=="  p__Crenarchaeota"| variable=="  p__Chloroflexi"| variable=="  p__Gemmatimonadetes"| variable=="  p__Planctomycetes"| variable=="  p__Firmicutes"| variable=="  p__Nitrospirae"| variable=="  c__Alphaproteobacteria"| variable=="  c__Betaproteobacteria"| variable=="  c__Deltaproteobacteria"| variable=="  c__Gammaproteobacteria")

head(dataset_comb_melt)
dataset_comb_melt$value<-dataset_comb_melt$value/10903

stats<-ddply(dataset_comb_melt, .(SoilFrac,variable),summarise,
means=mean(value),
se=sd(value)/sqrt(length(value))
)
final_int<-rbind(final_int, stats)

print(t/100)
}

final_int_summ<-ddply(subset(final_int, variable != "  p__"), .(SoilFrac,variable),summarise,
abundance=mean(means),
se_mean=mean(se)
)
library(ggplot2)

final_int_summ$variable<-factor(final_int_summ$variable,levels=rev(c("  p__Acidobacteria","  p__Verrucomicrobia","  p__Actinobacteria","  p__Bacteroidetes","  p__Planctomycetes","  p__Crenarchaeota","  p__Gemmatimonadetes","  p__Chloroflexi","  p__Firmicutes","  p__Nitrospirae","  c__Alphaproteobacteria","  c__Betaproteobacteria","  c__Deltaproteobacteria","  c__Gammaproteobacteria")))
levels(final_int_summ$variable)<-c("Gamma-Proteobacteria","Delta-Proteobacteria","Beta-Proteobacteria","Alpha-Proteobacteria","Nitrospirae","Firmicutes","Chloroflexi","Gemmatimonadetes","Crenarcheota","Planctomycetes","Bacteroidetes","Actinobacteria","Verrucomicrobia","Acidobacteria")
final_int_summ$SoilFrac<-factor(final_int_summ$SoilFrac)

final_int_summ$SoilFrac<-factor(final_int_summ$SoilFrac,levels=rev(c("Micro","SM","MM","LM","WS")))
levels(final_int_summ$SoilFrac)
levels(final_int_summ$SoilFrac)<-rev(c("Micro","Small","Medium","Large","Whole Soil"))
ggplot(final_int_summ)+geom_pointrange(aes(x=variable,y=abundance*100,ymax=(abundance+se_mean)*100,ymin=(abundance-se_mean)*100,colour=SoilFrac),position=position_dodge(width=.9),shape=32)+geom_bar(aes(x=variable, y=abundance*100,fill=SoilFrac),stat="identity",position="dodge")+coord_flip()+theme_bw(base_size=15)+theme(aspect.ratio=1)+scale_y_log10()+scale_colour_manual(name="Soil\nFraction",values=rev(c("black","#d6604d","#b2182b","#67001f","#878787")))+scale_fill_manual(name="Soil\nFraction",values=rev(c("black","#d6604d","#b2182b","#67001f","#878787")))+labs(y="Relative Abundance (%)",x="Taxonomic Group")


# looking at the statistics
stat_results<-read.table("~/Google Drive/postdoc_writing/COBS_16S_round3/phylum_statistic_results.txt",header=T,sep=" ")

head(stat_results)
ddply(subset(stat_results,trt=="SoilFrac"),.(taxon),summarise,
mean_p=mean(Pr..F.),
hi=quantile(Pr..F., 0.975)
)

summary(lm(c(1,2,3,4)~c(5,6,7,8)))