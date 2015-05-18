frac_div<-read.table("~/Google Drive/postdoc_writing/COBS_16S_round3/fracs_diversity_estimate_results.txt",header=T,sep=" ")

head(frac_div)

stats<-ddply(frac_div, .(Crop,SoilFrac,variable),summarise,
value=mean(val),
SE=mean(se)
)


stats$SoilFrac<-factor(stats$SoilFrac, levels=c("Micro","SM","MM","LM","WS"))
levels(stats$SoilFrac)<-c("Micro","Small","Medium","Large","Whole\nSoil")
stats$variable<-factor(stats$variable, levels=c("alpha","H","E"))
levels(stats$variable)<-c("Richness","Shannon's Diversity (H')","Evenness")

ggplot(stats)+geom_pointrange(aes(x=SoilFrac, y=value, ymax=value+SE,ymin=value-SE,colour=SoilFrac, shape=Crop),size=.6,position=position_dodge(width=0.7))+geom_line(aes(x=SoilFrac,y=value,group=Crop,linetype=Crop),alpha=0.5)+facet_wrap(~variable, scales="free")+theme_bw(base_size=13)+theme(aspect.ratio=1,axis.text.x=element_text(angle=25))+scale_colour_manual(values=c("black","#d6604d","#b2182b","#67001f","#878787"))+guides(colour=F)+scale_shape_discrete(name="Cropping\nSystem",labels=c("Corn","Prairie","Fertilized\nPrairie"))+labs(x="Soil Fraction", y="Mean Value")+scale_linetype_manual(name="Cropping\nSystem",labels=c("Corn","Prairie","Fertilized\nPrairie"),values=c("solid","dashed","dotted"))
