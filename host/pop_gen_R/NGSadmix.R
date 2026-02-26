setwd('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/pop_gen/NGSadmix/')

##### read in data #####
# sftped qopt and log files from HPC
# ran 20 iterations of each K (1-7)
# adapted from: https://baylab.github.io/MarineGenomics/week-9-population-structure-using-ngsadmix.html

# set directory where log and qopt files are
dir="bams90"
# set prefix for log and qopt files
prefix="bams90_k"
title=dir #saving title for plotting

#read in the data
data<-list.files(dir,pattern = ".log", full.names = T)

#use lapply to read in all our log files at once
# change 1:XX to number of iterations*Ks
bigData<-lapply(1:140, FUN = function(i) readLines(data[i]))
# find the line that starts with "best like=" or just "b"
library(stringr)
#this will pull out the line that starts with "b" from each file and return it as a list
foundset<-sapply(1:140, FUN= function(x) bigData[[x]][which(str_sub(bigData[[x]], 1, 1) == 'b')])
#now we need to pull out the first number in the string (max likelihood value), we'll do this with the function sub
as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )
#now lets store it in a dataframe

#make a dataframe with an index 1:10 for 20 iterations, this corresponds to our K values
logs<-data.frame(K = rep(1:7, each=20))
#add to it our likelihood values
logs$like<-sapply(strsplit(split='=',as.character(foundset)), FUN='[',2)
logs$like<-as.numeric(sapply(strsplit(split=' after ',logs$like), FUN='[',1))

##### choosing a K ####
## get an idea of variation across runs, the lower the better
tapply(logs$like, logs$K, FUN= function(x) sd(abs(x))/mean(abs(x)))

## following Evanno 2005 https://onlinelibrary.wiley.com/doi/10.1111/j.1365-294X.2005.02553.x
library(ggplot2)
library(dplyr)
likes=as.data.frame(logs %>% group_by(K) %>% summarize(mean=mean(like),sd=sd(like)))
plot1=ggplot(likes,aes(x=K,y=mean))+geom_point()+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=1)+ylab('L(k)')+
  geom_line()
# rate of change (Lprime) is calculated as L(k)-L(k-1) aka the difference in likelihoods between Ks
Lprime=c()
K=c(1:7)
for (i in K) {
  delL=mean(logs[logs$K==i,]$like) - mean(logs[logs$K==(i-1),]$like)
  Lprime=c(Lprime,delL)
}
Lprime=as.data.frame(cbind(K,Lprime))
plot2=ggplot(Lprime,aes(x=K,y=Lprime))+geom_point()+geom_line()+ylab('L\'(k)')

# now second order rate of change (Lpp) which is Lprime(L+1)-Lprime(k)
Lpp=c()
K=c(1:7)
for (i in K) {
  delLprime=abs(Lprime[Lprime$K==(i+1),]$Lprime - Lprime[Lprime$K==(i),]$Lprime)
  Lpp=c(Lpp,delLprime)
}
Lpp=as.data.frame(cbind(K,Lpp))
plot3=ggplot(Lpp,aes(x=K,y=Lpp))+geom_point()+geom_line()+ylab('L\'\'(k)')
# Nas OK here, cannot estimate for K=1 or highest K (becasue no K+1)

# finally delta k which is Lpp / sd(like)
sd=tapply(logs$like, logs$K, FUN= function(x) sd(abs(x)))
delK=as.data.frame(cbind(Lpp,sd))
delK$delK=abs(delK$Lpp)/abs(delK$sd)
plot4=ggplot(delK,aes(x=K,y=delK))+geom_point()+geom_line()+ylab('delK')

library(gridExtra)
library(grid)
grid.arrange(plot1,plot2,plot3,plot4,ncol=2,top=textGrob(title))

##### Plotting Ks ######
plot_list = list()
inds=read.table("bams90/bams90") # change to list of individuals
i2p=cbind(inds,sapply(strsplit(inds$V1,split='_'),FUN='[',1))
names(i2p)=c("ind","pop")

##### change pop order to match order you want sites to be (and number if subsetting)
pop_order=c('OK','MI','TK','LK','BP','SI')
#pop_order=c('OK','MI','TK','LK','BP')
#pop_order=c('SI')

for (K in 1:6){
Kdf=read.table(paste(dir,"/",prefix,K,"run10.qopt",sep=""),header=F)

Kdf=cbind(Kdf,i2p)
row.names(Kdf)=Kdf$ind

library(tidyr)

Kdf$pop=factor(Kdf$pop,levels=pop_order)
Kdf_order=Kdf[order(Kdf$pop,Kdf$V1),]
sample_order=Kdf_order$ind

K_melt <- Kdf %>% gather(key=admix,value=frac,-pop,-ind)
K_melt$pop=factor(K_melt$pop,levels=pop_order)
K_melt$ind=factor(K_melt$ind,levels=sample_order)
K_order=K_melt[order(K_melt$ind),]
labs=pop_order

library(forcats)
plot=ggplot(K_order, aes(fill=admix, y=frac, x=ind)) + 
  geom_bar(position="stack", stat="identity",width=1) +
  ylab('ancestry')+
  xlab(NULL)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
  theme_bw(base_size = 18) + theme(legend.position = "none") + theme(axis.text.x=element_blank()) +
  facet_grid(~fct_inorder(pop), scales = "free", switch = "x", space = "free") +
  theme(panel.spacing.x = grid:::unit(0, "lines"), 
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.ticks.x = element_blank())+
  ggtitle(paste("K =",K))

plot_list[[K]] = plot
}

#stitch together
grid.arrange(plot_list[[2]]+theme(plot.margin=unit(c(0.1,0.5,0.1,0),"cm"),plot.title=element_text(size=16),strip.text.x = element_blank())+scale_y_continuous(breaks=seq(0,1,by=0.5),expand = c(0, 0))+ylab(''),
             plot_list[[3]]+theme(plot.margin=unit(c(0.1,0.5,0.1,0),"cm"),plot.title=element_text(size=16),strip.text.x = element_blank())+scale_y_continuous(breaks=seq(0,1,by=0.5),expand = c(0, 0))+ylab(''),
             plot_list[[4]]+theme(plot.margin=unit(c(0.1,0.5,0.1,0),"cm"),plot.title=element_text(size=16),strip.text.x = element_blank())+scale_y_continuous(breaks=seq(0,1,by=0.5),expand = c(0, 0))+ylab(''),
             plot_list[[5]]+theme(plot.margin=unit(c(0.1,0.5,0.1,0),"cm"),plot.title=element_text(size=16),strip.text.x = element_blank())+scale_y_continuous(breaks=seq(0,1,by=0.5),expand = c(0, 0))+ylab(''),
             plot_list[[6]]+theme(plot.margin=unit(c(0.1,0.5,0.1,0),"cm"),plot.title=element_text(size=16),strip.text.x = element_blank())+scale_y_continuous(breaks=seq(0,1,by=0.5),expand = c(0, 0))+ylab(''),
             ncol=1)

###### Plot Best K #######
K=2 # set to best K
Kdf=read.table(paste(dir,"/",prefix,K,"run15.qopt",sep=""),header=F)

Kdf=cbind(Kdf,i2p)
row.names(Kdf)=Kdf$ind
Kdf$pop=factor(Kdf$pop,levels=pop_order)
Kdf_order=Kdf[order(Kdf$pop,Kdf$V1,Kdf$V2),]
sample_order=Kdf_order$ind

K_melt <- Kdf %>% gather(key=admix,value=frac,-pop,-ind)
K_melt$pop=factor(K_melt$pop,levels=pop_order)
K_melt$ind=factor(K_melt$ind,levels=sample_order)
K_order=K_melt[order(K_melt$ind),]
saveRDS(K_order,file=paste('final/dfs/',prefix,K,'.Rdata',sep=''))

labs=pop_order

ggplot(K_order, aes(fill=admix, y=frac, x=ind)) + 
  geom_bar(position="stack", stat="identity",width=1) +
  ylab('ancestry')+
  xlab(NULL)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
  theme_bw(base_size = 18) + theme(legend.position = "none") + theme(axis.text.x=element_blank()) +
  facet_grid(~fct_inorder(pop), scales = "free", switch = "x", space = "free") +
  theme(panel.spacing.x = grid:::unit(0, "lines"), 
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.ticks.x = element_blank())

#### how many admixed ####
# for all genets
# in this case, SI is assigned mostly to V1
# check how many inds from other pops have >25% assignment to V1
Kdf %>% group_by(pop) %>% summarise(sum(V1>0.25))






