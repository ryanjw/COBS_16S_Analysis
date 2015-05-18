ws1<-read.table("~/Google Drive/postdoc_writing/COBS_16S_round3/WS_cocur_results/WS_cocur_results.txt", header=T,sep=" ")
ws2<-read.table("~/Google Drive/postdoc_writing/COBS_16S_round3/WS_cocur_results/WS_cocur_results2.txt", header=T,sep=" ")
ws3<-read.table("~/Google Drive/postdoc_writing/COBS_16S_round3/WS_cocur_results/WS_cocur_results3.txt", header=T,sep=" ")
ws4<-read.table("~/Google Drive/postdoc_writing/COBS_16S_round3/WS_cocur_results/WS_cocur_results4.txt", header=T,sep=" ")

names(ws1)[1:4]<-c("a","SoilFrac","n1","n2")
names(ws2)[1:4]<-c("a","SoilFrac","n1","n2")
names(ws3)[1:4]<-c("a","SoilFrac","n1","n2")
names(ws4)[1:4]<-c("a","SoilFrac","n1","n2")

library(fdrtool)

as<-as.vector(unique(ws1$a))
ws1q<-data.frame()
for(i in 1:length(as)){
	temp<-subset(ws1, a==as[i])
	temp$qval<-fdrtool(temp$p.value,statistic="pvalue",plot=F,verbose=F)$qval
	ws1q<-rbind(ws1q,temp)
}
as<-as.vector(unique(ws2$a))
ws2q<-data.frame()
for(i in 1:length(as)){
	temp<-subset(ws2, a==as[i])
	temp$qval<-fdrtool(temp$p.value,statistic="pvalue",plot=F,verbose=F)$qval
	ws2q<-rbind(ws2q,temp)
}
as<-as.vector(unique(ws3$a))
ws3q<-data.frame()
for(i in 1:length(as)){
	temp<-subset(ws3, a==as[i])
	temp$qval<-fdrtool(temp$p.value,statistic="pvalue",plot=F,verbose=F)$qval
	ws3q<-rbind(ws3q,temp)
}
as<-as.vector(unique(ws4$a))
ws4q<-data.frame()
for(i in 1:length(as)){
	temp<-subset(ws4, a==as[i])
	temp$qval<-fdrtool(temp$p.value,statistic="pvalue",plot=F,verbose=F)$qval
	ws4q<-rbind(ws4q,temp)
}

ws_all<-rbind(ws1q,ws2q,ws3q,ws4q)

head(ws_all)

dataset<-subset(ws_all,qval < 0.01 & n1 != "  g__" & n2 != "  g__" & rho > 0)
library(plyr)
dataset$pair<-paste(dataset$n1,dataset$n2)
data_av2<-ddply(dataset, .(SoilFrac,n1,n2,pair),summarise,.progress="text",
av_rho=mean(rho),
hi=quantile(rho,0.975),
lo=quantile(rho,0.025),
SD=sd(rho)
)
library(igraph)
fracs<-as.vector(unique(data_av2$SoilFrac))
dg_final<-data.frame()
btwn_final<-data.frame()
ratios<-data.frame()
edge_btwn_final<-data.frame()

# loop for producing network stats
for(i in 1:length(fracs)){
	#i<-1
	temp<-subset(data_av2, SoilFrac==fracs[i])
	g_temp<-graph.edgelist(as.matrix(temp[,2:3]),directed=FALSE)
	E(g_temp)$weight<-temp$av_rho
	g_temp<-simplify(g_temp)
	g_temp<-delete.vertices(g_temp,which(degree(g_temp) < 1))
	
	#clustering ratio
	clustering_coeff<-transitivity(g_temp,type="global")
	clustering_coeff_rand<-transitivity(erdos.renyi.game(length(V(g_temp)),length(E(g_temp)),type="gnm"))
	cluster_ratio<-clustering_coeff/clustering_coeff_rand
	ratio<-data.frame(fracs[i],cluster_ratio)
	names(ratio)[1]<-"SoilFrac"
	ratios<-rbind(ratios,ratio)
	
	#degree
	dg_df<-data.frame(degree(g_temp, normalized=T))
	dg_df<-cbind(rownames(dg_df),dg_df)
	rownames(dg_df)<-NULL
	names(dg_df)<-c("g","norm_degree")
	dg_df<-arrange(dg_df, -norm_degree)
	dg_df<-dg_df[1:5,]
	dg_df$SoilFrac<-fracs[i]
	dg_final<-rbind(dg_final, dg_df)
	
	#betweenness
	btwn_df<-data.frame(betweenness(g_temp, normalized=T,weights=E(g_temp)$weight,directed=F))
	btwn_df<-cbind(rownames(btwn_df),btwn_df)
	rownames(btwn_df)<-NULL
	names(btwn_df)<-c("g","norm_btwn")
	btwn_df<-arrange(btwn_df,-norm_btwn)
	btwn_df<-btwn_df[1:5,]
	btwn_df$SoilFrac<-fracs[i]
	btwn_final<-rbind(btwn_final,btwn_df)
	
	# edge_betweenness
	edge_btwn<-data.frame(get.edgelist(g_temp),edge.betweenness(g_temp,e=E(g_temp),directed=F,weights=E(g_temp)$weight))
	names(edge_btwn)<-c("g1","g2","edge_btwn")
	edge_btwn$SoilFrac<-fracs[i]
	edge_btwn<-arrange(edge_btwn, -edge_btwn)
	edge_btwn<-edge_btwn[1:5,]
	edge_btwn_final<-rbind(edge_btwn_final,edge_btwn)
		
}
dg_final
btwn_final
ratios

edge_btwn_final

# plotting
library(GGally)
library(intergraph)
head(data_av2)
g<-graph.edgelist(as.matrix(subset(data_av2)[,2:3]),directed=FALSE)
E(g)$weight<-subset(data_av2)$av_rho
g<-simplify(g)
g<-delete.vertices(g,which(degree(g) < 1))
otus<-read.csv("~/Google Drive/COBS_16S_biodiversity/otu_list.csv",header=T,check.names=F)
otus<-subset(otus, genus != "  g__")
genera_list<-data.frame()
generas<-as.vector(unique(otus$genus))
for(i in 1:length(generas)){
	#i<-1
	temp<-subset(otus, genus==generas[i])
	temp<-temp[1,-c(1,8)]
	genera_list<-rbind(genera_list,temp)
	print(i/length(generas))
}

chems<-data.frame(rbind(rep("Extractable.C.gC.m.2",6),rep("Extractable.N.gC.m.2",6)))
names(chems)<-names(genera_list)
genera_list<-rbind(genera_list,chems)
gnetwork<-asNetwork(g)
thing<-asDF(gnetwork)
df_g<-thing$vertexes
dim(df_g)
head(genera_list)
df_g2<-data.frame()
for(i in 1:dim(df_g)[1]){
	r<-df_g[i,]
	list_row<-subset(genera_list, genus==r$vertex.names[1] )[1,]
	r<-merge(r,list_row,by.x="vertex.names",by.y="genus")
	df_g2<-rbind(df_g2,r)
}
head(df_g2)
dim(df_g2)
df_g<-df_g2
df_g<-arrange(df_g,intergraph_id)
df_g$categories<-"Other"
for(i in 1:dim(df_g)[1]){
	if(df_g$phylum[i]=="  p__Actinobacteria"){df_g$categories[i]<-"Actinobacteria"}
	if(df_g$class[i]=="  c__Alphaproteobacteria"){df_g$categories[i]<-"Alphaproteobacteria"}
	if(df_g$class[i]=="  c__Betaproteobacteria"){df_g$categories[i]<-"Betaproteobacteria"}
	if(df_g$phylum[i]=="  p__Firmicutes"){df_g$categories[i]<-"Firmicutes"}
	if(df_g$class[i]=="  c__Gammaproteobacteria"){df_g$categories[i]<-"Gammaproteobacteria"}
	if(df_g$phylum[i]=="  p__Bacteroidetes"){df_g$categories[i]<-"Bacteroidetes"}
	if(df_g$class[i]=="  c__Deltaproteobacteria"){df_g$categories[i]<-"Deltaproteobacteria"}
	if(df_g$phylum[i]=="  p__Verrucomicrobia"){df_g$categories[i]<-"Verrucomicrobia"}
	if(df_g$phylum[i]=="Extractable.C.gC.m.2"){df_g$categories[i]<-"Total Carbon"}	
	if(df_g$phylum[i]=="Extractable.N.gC.m.2"){df_g$categories[i]<-"Total Nitrogen"}		
}
library(ggplot2)
library(GGally)
df_g
ws_network<-ggnet(g,size=0,method="kamadakawaii",segment.size=log(1+thing$edges$weight))+geom_point(aes(colour=df_g$categories,shape=df_g$categories),size=5)+scale_colour_manual(name="Taxonomic\nGroups",values=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","blue4","#999999"))+scale_shape_manual(name="Taxonomic\nGroups",values=c(rep(16,8),rep(17,1),16))+theme(text=element_text(size=17))+guides(color=F,shape=F)


