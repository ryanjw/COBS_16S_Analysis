hist(rpois(100,10000))
hist(round(rgamma(1000,1,rate=1))*1000)

# simulating communities
fseries<-function(X,C){
	y<-(-1/log10(1-C))*(C^X/X)
	return(y)
}

rank<-seq(1,10000,1)
df<-data.frame(rank)
head(df)
df$ab<-fseries(df$rank, .96)
df$rel_ab<-df$ab/sum(df$ab)

dim(subset(df, ab > 0))

cs<-seq(.9,1,.0001)
df<-data.frame()
for(i in 1:length(cs)){
	#i<-1
	fs<-(fseries(rank,cs[i]))
	alpha<-length(fs[fs>0])
	dftemp<-data.frame(cs[i],alpha)
	df<-rbind(df,dftemp)
}
ggplot(df)+geom_point(aes(x=cs.i.,y=alpha))+geom_hline(yintercept=10000)
head(subset(df, alpha >= 10000))

# for this many otus, need to start at .9291
C<-seq(.9300,.9999,length.out=100)
C

#How does C affect evenness and richness
shannon<-function(XXX){
	denom<-1
	for(elem in 1:length(XXX)){
		denom<-denom*XXX[elem]^XXX[elem]
	}
	return(log(1/denom))
}

df_final<-data.frame()
for(i in 1:length(C)){
	#i<-1
	df<-data.frame(rank)
	df$ab<-(fseries(rank,C[i]))
	df$rel_ab<-df$ab/sum(df$ab)
	alpha<-specnumber(df$rel_ab,MARGIN=2)
	shannon<-diversity(df$rel_ab,MARGIN=2)
	evenness<-shannon/log(alpha)
	new_row<-data.frame(C[i],alpha,shannon,evenness)
	df_final<-rbind(df_final,new_row)
}

head(df_final)
library(reshape2)
df_final_melt<-melt(df_final, id=c("C.i."))
ggplot(df_final_melt)+geom_point(aes(x=C.i.,y=value))+facet_wrap(~variable,scales="free")

#shannons and evenness changes

reads<-seq(5000,20000,length.out=100)
C<-seq(.9300,.9999,length.out=100)
df_final<-data.frame()
for(r in 1:length(reads)){
	#r<-1
	for(i in 1:length(C)){
		#i<-1
		df<-data.frame(rank)
		df$ab<-(fseries(rank,C[i]))
		df$rel_ab<-df$ab/sum(df$ab)
		collection<-sample(df$rank,reads[r],replace=T,prob=df$rel_ab)
		community<-count(collection)
		alpha<-specnumber(community$x)
		new_row<-data.frame(reads[r],C[i],alpha)
		df_final<-rbind(df_final,new_row)
	}
	print(r/length(reads))
}
names(df_final)<-c("reads","C","alpha")
ggplot(df_final)+geom_tile(aes(x=reads,y=C,fill=(alpha)))+scale_fill_gradient(low="red",high="darkblue")



# first generate species list

# let's generate a matrix with 10,000 total unique OTUs
otus<-matrix(ncol=1,nrow=0)
for(i in 1:10000){
	otus<-rbind(otus, paste(sample(letters,6),collapse=""))
}

rank<-seq(1,10000,1)
df<-data.frame(otus,rank)
head(df)
df$rel_ab<-fseries(df$rank,.99)
df$rel_ab<-df$rel_ab/sum(df$rel_ab)
