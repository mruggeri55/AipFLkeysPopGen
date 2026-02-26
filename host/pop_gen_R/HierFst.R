library(adegenet)
library(pegas)
library(vcfR)
library(hierfstat)
library(poppr)
setwd('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/pop_gen/')

vcf=read.vcfR('MLLs89/MLLs89.unlinked.bcf.gz')
bams=read.table("MLLs89/MLLs89",header=F)[,1]
bams=sub("\\.bam","",bams,perl=T)

genind <- vcfR2genind(vcf)
bams=as.data.frame(bams)
colnames(bams)='INDIVIDUAL'
bams$STRATA=sapply(strsplit(x=bams$INDIVIDUAL,split='_'),FUN='[',1)

strata(genind) <- bams
setPop(genind) <- ~STRATA

Fst=pairwise.WCfst(genind,diploid=TRUE)
Fst_order=Fst[c('OK','MI','TK','LK','BP','SI'),c('OK','MI','TK','LK','BP','SI'),drop = FALSE]
diag(Fst_order)<-0
write.csv(Fst_order,"/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/pop_gen/MLLs89/Fst_mat.csv")

library(corrplot)
quartz()
corrplot(Fst_order,type='lower', is.corr = F,diag = F,addCoef.col = 'black',
         tl.col = 'black', tl.cex=1.5,
         col.lim=c(0,0.1))
quartz.save('MLLs89/plots/HierFstat_Fst_corrplot_MLLs89.pdf',type='pdf')
dev.off()

## IBD
WaterDist=read.csv('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/gps/over_water_dist.csv')
WaterDist$PopPair=paste(WaterDist$Site1,WaterDist$Site2, sep = '_')
library(dplyr)
Dist = WaterDist %>% group_by(PopPair) %>% summarise(waterdist_km=mean(over_water_dist_km))
Dist=as.data.frame(Dist[!Dist$waterdist_km==0,])

head(Dist)
Dist$pop1=sapply(strsplit(Dist$PopPair,split='_'),FUN='[',1)
Dist$pop2=sapply(strsplit(Dist$PopPair,split='_'),FUN='[',2)
library(reshape2)
Distma=dcast(Dist,pop1~pop2,value.var = 'waterdist_km')
rownames(Distma)=Distma$pop1
Distma$pop1=NULL
Distma_order=Distma[c('OK','MI','TK','LK','BP','SI'),c('OK','MI','TK','LK','BP','SI'),drop = FALSE]
Distma=as.matrix(Distma_order)

plot(Distma,Fst_order)
melt1=reshape2::melt(Distma,dimnames(Distma),varnames = c('pop1','pop2'),value.name = 'waterdist_km')
melt2=reshape2::melt(Fst_order,dimnames(Fst_order),varnames = c('pop1','pop2'),value.name = 'Fst')
HierFst_long=merge(melt1,melt2,by=c('pop1','pop2'))
HierFst_long$Fst_norm=HierFst_long$Fst/(1 - HierFst_long$Fst)

# all
library(ggplot2)
ggplot(HierFst_long[!duplicated(HierFst_long$waterdist_km),],
       aes(x=waterdist_km,y=Fst_norm,label=paste(pop1,pop2)))+
  geom_point()+geom_text()+geom_smooth(method = "lm", se = FALSE) 
# remove OK comps
ggplot(HierFst_long[!duplicated(HierFst_long$waterdist_km) & !HierFst_long$pop1=='OK' & !HierFst_long$pop2=='OK',],
       aes(x=waterdist_km,y=Fst_norm,label=paste(pop1,pop2)))+
  geom_point()+geom_text()+geom_smooth(method = "lm", se = FALSE) 
# remove SIs
ggplot(HierFst_long[!duplicated(HierFst_long$waterdist_km) & !HierFst_long$pop1=='SI' & !HierFst_long$pop2=='SI',],
       aes(x=waterdist_km,y=Fst_norm,label=paste(pop1,pop2)))+
  geom_point()+geom_text()+geom_smooth(method = "lm", se = FALSE) 
# remove SI and OK
ggplot(HierFst_long[!duplicated(HierFst_long$waterdist_km) & !HierFst_long$pop1 %in% c('SI','OK') & !HierFst_long$pop2 %in% c('SI','OK'),],
       aes(x=waterdist_km,y=Fst_norm,label=paste(pop1,pop2)))+
  geom_point()+geom_text()+geom_smooth(method = "lm", se = FALSE) 
# OK only 
ggplot(HierFst_long[HierFst_long$pop1=='OK',],
       aes(x=waterdist_km,y=Fst_norm,label=pop2))+
  geom_point()+geom_text()+geom_smooth(method = "lm", se = FALSE)
#SI only
ggplot(HierFst_long[HierFst_long$pop1=='SI',],
       aes(x=waterdist_km,y=Fst_norm,label=pop2))+
  geom_point()+geom_text()+geom_smooth(method = "lm", se = FALSE)

# plot all together with different regression lines
HierFst_long$group=rep('G1',times=length(HierFst_long$pop1))
HierFst_long$group[HierFst_long$pop1 %in% c('MI','TK','LK','BP')]='G2'
HierFst_long$group[HierFst_long$pop1=='SI']='G3'
HierFst_long$group[HierFst_long$pop2=='SI']='G3'
HierFst_long$group[HierFst_long$pop1=='OK']='G1'
HierFst_long$group[HierFst_long$pop2=='OK']='G1'

ggplot(HierFst_long[!duplicated(HierFst_long$waterdist_km),],
       aes(x=waterdist_km,y=Fst_norm,label=paste(pop1,pop2,sep='-'),group=group))+
  geom_point()+
  geom_smooth(method = "lm", se = FALSE,aes(color=group))+
  geom_text(vjust=-0.25)+
  theme_bw(base_size=20)+
  theme(legend.position = 'none')+
  ylab("Fst/1-Fst")+
  xlab('over water distance (km)')+
  xlim(25,475)+
  ylim(0.01,0.08)

#### mantel tests ####
IBDma=dcast(HierFst_long, pop1 ~ pop2, value.var ='Fst_norm')
rownames(IBDma)=IBDma[,1]
IBDma[,1]=NULL

mantel.rtest(dist(IBDma),dist(Distma),nrepet=999)
#p=0.074

# no SIs
mantel.rtest(dist(IBDma[!rownames(IBDma)=='SI',!colnames(IBDma)=='SI']),
             dist(Distma[!rownames(Distma)=='SI',!colnames(Distma)=='SI']),nrepet=999)
#p=0.005

# no OKs
mantel.rtest(dist(IBDma[!rownames(IBDma)=='OK',!colnames(IBDma)=='OK']),
             dist(Distma[!rownames(Distma)=='OK',!colnames(Distma)=='OK']),nrepet=999)
# p=0.201 so not significant without OKs

# no OKs or SIs
mantel.rtest(dist(IBDma[!rownames(IBDma) %in% c('OK','SI'),!colnames(IBDma) %in% c('OK','SI')]),
             dist(Distma[!rownames(Distma) %in% c('OK','SI'),!colnames(Distma) %in% c('OK','SI')]),nrepet=999)
# p=0.111

# SIs only
mantel.rtest(dist(IBDma[rownames(IBDma)=='SI',!colnames(IBDma)=='SI']),
             dist(Distma[!rownames(IBDma)=='SI',colnames(Distma)=='SI']),nrepet=999)
# this doesnt work because not a pairwise matrix, try linear model
summary(lm(IBDma$SI[!colnames(IBDma)=='OK']~Distma_order$SI[!colnames(Distma_order)=='OK']))
#p=0.02, R2=0.9375, slope= -1.618e-04, intercept= 7.793e-02

# OK only
summary(lm(IBDma$OK~Distma_order$OK))
#p=0.001356, R2=0.9712, slope=-2.78e-04, intercept=1.579e-01


basicstat <- basic.stats(genind, diploid = TRUE, digits = 2) 
names(basicstat)  
basicstat$overall

# can bootstrap Fis
boot.ppfis(genind) 

Fis = as.data.frame(basicstat$Fis)
Ho = as.data.frame(basicstat$Ho)
Hs = as.data.frame(basicstat$Hs)

library(tidyr)
Fis_melt=Fis %>% gather(key=pop,value=Fis)
Ho_melt=Ho %>% gather(key=pop,value=Ho)

ggplot(Fis_melt,aes(x=pop,y=Fis))+geom_boxplot()
ggplot(Ho_melt,aes(x=pop,y=Ho))+geom_boxplot()
