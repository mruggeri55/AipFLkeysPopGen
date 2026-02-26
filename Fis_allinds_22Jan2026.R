setwd("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/vcfs")

library(vcfR)
library(adegenet)
vcf=read.vcfR('bams300.bcf.gz',convertNA = TRUE) 
genind = vcfR2genind(vcf)
genclone=as.genclone(genind)
genclone
nmll(genclone)

# add metadata
strata=as.data.frame(cbind(IND=labels(prev_dist),POP=sapply(strsplit(labels(prev_dist),'_'),FUN='[',1),
                           SUB=sapply(strsplit(labels(prev_dist),'_rt'),FUN='[',1),
                           RT=paste(sapply(strsplit(labels(prev_dist),'_rt'),FUN='[',1),sapply(strsplit(labels(prev_dist),'_'),FUN='[',3),sep='_')),
                     col.names=c('IND','POP','SUB','RT'))
strata(genclone) <- strata
setPop(genclone) <- ~POP

keep=indNames(genclone)[!indNames(genclone) %in% c(TR2s,'CC7_host.bam')]
gensub=genclone[keep,]

library(hierfstat)
basicstat <- basic.stats(gensub, diploid = TRUE, digits = 2) 
basicstat$overall
Ho <- colMeans(basicstat$Ho,na.rm=T)
He <- colMeans(basicstat$Hs,na.rm=T)
Fis<-colMeans(basicstat$Fis,na.rm = T)
gen_stats <- as.data.frame(cbind(Ho,He,Fis))
gen_stats$Site=row.names(gen_stats)

Fis = as.data.frame(basicstat$Fis)
library(reshape2)
Fis_long=melt(Fis,value.name='Fis',variable.name = 'pop')
Fis_long$pop=factor(Fis_long$pop,levels=c('OK','MI','TK','LK','BP','SI'))

library(ggplot2)
library(ggridges)
site_cols=readRDS('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/MS/final_figs/colPalSites.RDS')

ggplot(Fis_long,aes(x=pop,y=Fis,fill=pop))+
  geom_boxplot()+
  scale_fill_manual(values=site_cols$site)

ggplot(Fis_long,aes(x=pop,y=Fis,fill=pop))+
  geom_violin()+
  geom_boxplot(width=0.1)+
  scale_fill_manual(values=site_cols$site)+
  geom_hline(yintercept=0,linetype="dashed")+
  theme_bw(base_size=20)+
  labs(x="sampling location (N \u2192 S)",y=expression("F"[IS]))+
  theme(legend.position='none')

ggplot(Fis_long,aes(x=Fis,y=-0.5))+
  geom_boxplot(aes(fill = pop),width=0.5) +
  geom_density(aes(x = Fis,fill=pop), inherit.aes = FALSE) +
  stat_boxplot(geom = "vline", aes(xintercept = 0), linetype='dashed')+
  theme_bw(base_size=20)+
  facet_grid(pop~.)+
  ylab('density')+
  theme(legend.position = 'none')

# testing whether significant differences in variance
library(car)
leveneTest(Fis~pop,data=Fis_long) # significant
library(dplyr)
# to do post hoc test, calculate residuals then do an anova on residuals
med <- Fis_long %>% group_by(pop) %>% 
  mutate(dat.med = ifelse(Fis,median(Fis, na.rm=TRUE), ifelse(Fis==NA, NA)))
med$dat.med.res<-abs(med$Fis-med$dat.med)
# Then we run an ANOVA, and post-hoc if necessary:
levene.dat.aov<-aov(dat.med.res~pop,med)
summary(levene.dat.aov)
levenePW=as.data.frame(TukeyHSD(levene.dat.aov)$pop)
# everything significantly different from OK
# everything signifcantly different from SI
# TK-MI different
write.csv(levenePW,"/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/Fis/results/LevenePW_allinds_22Jan2026.csv")

library(corrplot)
library(reshape2)
levenePW$pop1=sapply(strsplit(rownames(levenePW),'-'),FUN='[',1)
levenePW$pop2=sapply(strsplit(rownames(levenePW),'-'),FUN='[',2)
levenePW$p.adj=levenePW$`p adj`

cor_flat=subset(levenePW,select=c(pop1,pop2,diff,p.adj))
corr = dcast(cor_flat,pop1 ~ pop2, value.var ='diff')
rownames(corr)=corr$pop1
corr$pop1=NULL
corr=as.matrix(corr)
corr=corr[c('SI','BP','LK','TK','MI'),c('BP','LK','TK','MI','OK')]
corr[lower.tri(corr)] <- t(corr[upper.tri(corr)])
corr=corr[c('MI','TK','LK','BP','SI'),c('OK','MI','TK','LK','BP')]

p.mat = dcast(cor_flat,pop1 ~ pop2, value.var ='p.adj')
rownames(p.mat)=p.mat$pop1
p.mat$pop1=NULL
p.mat=as.matrix(p.mat)
p.mat=p.mat[c('SI','BP','LK','TK','MI'),c('BP','LK','TK','MI','OK')]
p.mat[lower.tri(p.mat)] <- t(corr[upper.tri(p.mat)])
p.mat=p.mat[c('MI','TK','LK','BP','SI'),c('OK','MI','TK','LK','BP')]

#quartz()
corrplot(corr,type='lower', is.corr = F,
         tl.col = 'black', tl.cex=1.5, 
         cl.ratio = 0.4, cl.cex = 0.5,
         col = COL2("RdYlBu", 12),
         p.mat=p.mat,insig = 'label_sig',
         sig.level = c(.001, .01, .05),pch.cex=1)
#quartz.save('Fis_corrplot.pdf',type='pdf')
#dev.off()
