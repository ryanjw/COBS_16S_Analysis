sizes<-fseries(c(1,2,3,4),.2)
reads<-10000
read_prop<-round(sizes/sum(sizes)*reads)
rank<-seq(1,10000,1)
fs<-fseries(rank,.999)
rel_ab<-fs/sum(fs)
full_sample<-c(sample(rank,read_prop[1],replace=T,prob=rel_ab),
sample(rank,read_prop[2],replace=T,prob=rel_ab),
sample(rank,read_prop[3],replace=T,prob=rel_ab),
sample(rank,read_prop[4],replace=T,prob=rel_ab))
comm<-count(full_sample)
alpha<-specnumber(comm[,2],MARGIN=2)
shannon<-diversity(comm[,2],MARGIN=2)
evenness<-shannon/log(alpha)
return(c(alpha,shannon,evenness))


fseries<-function(X,C){
	y<-(-1/log10(1-C))*(C^X/X)
	return(y)
}



comm_simm_diversity<-function(samp_skew,read_count,F_skewness){
	sizes<-fseries(c(1,2,3,4),samp_skew)
	reads<-read_count
	read_prop<-round(sizes/sum(sizes)*reads)
	rank<-seq(1,10000,1)
	fs<-fseries(rank,F_skewness)
	rel_ab<-fs/sum(fs)
	full_sample<-c(sample(rank,read_prop[1],replace=T,prob=rel_ab),
	sample(rank,read_prop[2],replace=T,prob=rel_ab),
	sample(rank,read_prop[3],replace=T,prob=rel_ab),
	sample(rank,read_prop[4],replace=T,prob=rel_ab))
	comm<-count(full_sample)
	alpha<-specnumber(comm[,2],MARGIN=2)
	shannon<-diversity(comm[,2],MARGIN=2)
	evenness<-shannon/log(alpha)
	return(c(alpha,shannon,evenness))
}

comm_simm_diversity(.5,10000,.9999)
skews<-seq(.1,.99999,length.out=50)
C<-seq(.93,.9999,length.out=50)
Reads<-seq(5000,20000,length.out=50)


final<-data.frame()
for(s in 1:length(skews)){
	for(cs in 1:length(C)){
		for(r in 1:length(Reads)){
			new_row<-data.frame(skews[s],C[cs],Reads[r],comm_simm_diversity(skews[s],Reads[r],C[cs]))
			final<-rbind(final,new_row)
		}
	}
	print(s/50)
}
?count
final<-read.table("~/Desktop/final_sim.txt",sep=" ",header=T)

head(final)

names(final)<-c("skew","C","reads","div_value")
final$div_variable<-rep(c("alpha","shannon","evenness"))
library(reshape2)
final_melt<-melt(final, id=c("div_variable","div_value"))
head(final_melt)
head(final)

ggplot(subset(final, div_variable=="alpha"))+geom_tile(aes(x=skew,y=C,fill=div_value))+scale_fill_gradient(low="red",high="darkblue",name="Alpha\nDiversity")+theme_bw(base_size=17)+theme(aspect.ratio=1)+labs(x="Micro-habitat Skewnness",y="Community Skewness")
ggplot(subset(final, div_variable=="alpha"))+geom_tile(aes(x=reads,y=C,fill=div_value))+scale_fill_gradient(low="red",high="darkblue",name="Alpha\nDiversity")+theme_bw(base_size=17)+theme(aspect.ratio=1)+labs(x="Sampling Effort",y="Community Skewness")
ggplot(subset(final, div_variable=="alpha"))+geom_tile(aes(x=reads,y=skew,fill=div_value))+scale_fill_gradient(low="red",high="darkblue",name="Alpha\nDiversity")+theme_bw(base_size=17)+theme(aspect.ratio=1)+labs(x="Sampling Effort",y="Micro-habitat Skewnness")


# making 4 different communities that span differences in sharing and difference in otus


# let's generate a matrix with 10,000 total unique OTUs
otus<-matrix(ncol=1,nrow=0)
for(i in 1:10000){
	otus<-rbind(otus, paste(sample(letters,6),collapse=""))
}

samp_skew<-.5
reads<-10000
F_skewness<-.9999
shared<-6000

sep_comm_func<-function(samp_skew,reads,F_skewness,shared){
left<-(10000-shared)/4

core<-otus[1:shared]
sub1<-sample(c(core,otus[(shared+1):(shared+left)]))
sub2<-sample(c(core,otus[(shared+left+1):(shared+left*2)]))
sub3<-sample(c(core,otus[(shared+left*2+1):(shared+left*3)]))
sub4<-sample(c(core,otus[(shared+left*3+1):(shared+left*4)]))

sizes<-fseries(c(1,2,3,4),samp_skew)
read_prop<-round(sizes/sum(sizes)*reads)
fs1<-fseries(seq(1,length(sub1),1),F_skewness)
fs2<-fseries(seq(1,length(sub2),1),F_skewness)
fs3<-fseries(seq(1,length(sub3),1),F_skewness)
fs4<-fseries(seq(1,length(sub4),1),F_skewness)
rel_ab1<-fs1/sum(fs1)
rel_ab2<-fs2/sum(fs2)
rel_ab3<-fs3/sum(fs3)
rel_ab4<-fs4/sum(fs4)
full_sample<-c(sample(sub1,read_prop[1],replace=T,prob=rel_ab1),
	sample(sub2,read_prop[2],replace=T,prob=rel_ab2),
	sample(sub3,read_prop[3],replace=T,prob=rel_ab3),
	sample(sub4,read_prop[4],replace=T,prob=rel_ab4))
full_sample_even<-c(sample(sub1,round(reads/4),replace=T,prob=rel_ab1),
	sample(sub2,round(reads/4),replace=T,prob=rel_ab2),
	sample(sub3,round(reads/4),replace=T,prob=rel_ab3),
	sample(sub4,round(reads/4),replace=T,prob=rel_ab4))
comm<-count(full_sample)
	alpha<-specnumber(comm[,2],MARGIN=2)
	shannon<-diversity(comm[,2],MARGIN=2)
	evenness<-shannon/log(alpha)
comm_even<-count(full_sample_even)
	alpha_even<-specnumber(comm_even[,2],MARGIN=2)
	shannon_even<-diversity(comm_even[,2],MARGIN=2)
	evenness_even<-shannon_even/log(alpha_even)
return(data.frame(samp_skew,reads,F_skewness,shared,alpha,shannon,evenness,alpha_even,shannon_even,evenness_even))

}

sep_comm_func(.5,10000,.9999,6000)

skews<-seq(.1,.99999,length.out=9)
C<-seq(.93,.9999,length.out=9)
Reads<-round(seq(5000,20000,length.out=9))
shareds<-seq(1000,9000,1000)

final<-data.frame()
for(s in 1:length(skews)){
	for(cs in 1:length(C)){
		for(r in 1:length(Reads)){
			for(k in 1:length(shareds)){
				for(i in 1:1000){
				new_row<-sep_comm_func(skews[s],Reads[r],C[cs],shareds[k])
				final<-rbind(final,new_row)}
				#print("Done")
			}
			
		}
	}
	print(s/9)
}