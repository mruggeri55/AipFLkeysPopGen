setwd("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/vcfs")

library(vcfR)
library(adegenet)
library(poppr)

vcf=read.vcfR('bams300.bcf.gz',convertNA = TRUE) 
genind = vcfR2genind(vcf)
genclone=as.genclone(genind)
genclone
nmll(genclone)

#bit_dist=bitwise.dist(genclone)
prev_dist=prevosti.dist(genclone)
ibs_dist=as.dist(read.table('../dist/bams300.ibsmat'))
#bruvo_dist=bruvo.dist(genclone)
hist(prev_dist)
#hist(bit_dist)
hist(ibs_dist)
#hist(bruvo_dist)

gac <- genotype_curve(genclone, sample = 100, maxloci=100, quiet = TRUE)

# quartz()
# genclone_filtered <- filter_stats(genclone, distance = prev_dist, plot = TRUE)
# #quartz.save('../poppr_strict/poppr2025_filterstats.pdf',type='pdf')
# dev.off()
genclone_filtered <- filter_stats(genclone, distance = prev_dist, plot = FALSE)


# In short, farthest neighbor is the most conservative, nearest neighbor can have a chaining effect, and average neighbor is somewhere in between. 
# Your choice of algorithm really depends on the biology of your organism.

# how to choose threshold?
# from poppr vignette:
# One method described in the literature of choosing a threshold is to look for an initial, small peak in the histogram of pairwise genetic distances and set the threshold to be between that peak and the larger peak `Bailleul, Stoeckel, and Arnaud-Haond (2016). This initial peak likely represents clones differentiated by a small set of random mutations.
# However, if this peak is not obvious, then another method is to look for the largest gap between all putative thresholds.

# our peak is pretty obvious
# let's see what the cutoff estimates say
print(farthest_thresh <- cutoff_predictor(genclone_filtered$farthest$THRESHOLDS))
print(average_thresh  <- cutoff_predictor(genclone_filtered$average$THRESHOLDS))
print(nearest_thresh  <- cutoff_predictor(genclone_filtered$nearest$THRESHOLDS))

# farthest neighbor is most conservative, let's use that
mlg.filter(genclone, distance = prev_dist, algorithm = "f") <- farthest_thresh
genclone
# add metadata
strata=as.data.frame(cbind(IND=labels(prev_dist),POP=sapply(strsplit(labels(prev_dist),'_'),FUN='[',1),
                           SUB=sapply(strsplit(labels(prev_dist),'_rt'),FUN='[',1),
                           RT=paste(sapply(strsplit(labels(prev_dist),'_rt'),FUN='[',1),sapply(strsplit(labels(prev_dist),'_'),FUN='[',3),sep='_')),
                           col.names=c('IND','POP','SUB','RT'))
strata(genclone) <- strata
setPop(genclone) <- ~POP

# technical rep distance
library(usedist)
TRs=labels(prev_dist)[grep("TR",labels(prev_dist))]
# excluding MI_sub3_rt4 because TR1 failed prep
# also adding in accidental tech reps
TRs=c(TRs[!grepl("MI_sub3_rt4",TRs)],'OK_sub2_rt3_A1_host.bam','OK_sub2_rt3_A2_host.bam','TK_sub1_rt3_93_host.bam','TK_subNA_rtNA_93B_host.bam')
subdist=as.matrix(dist_subset(prev_dist,TRs))
TR1s=c(TRs[grep('TR1',TRs)],'OK_sub2_rt3_A1_host.bam','TK_sub1_rt3_93_host.bam')
TR2s=c(TRs[grep('TR2',TRs)],'OK_sub2_rt3_A2_host.bam','TK_subNA_rtNA_93B_host.bam')
TRpw=subdist[row.names(subdist) %in% TR1s,colnames(subdist) %in% TR2s]
TRpw_dist=diag(TRpw)
names(TRpw_dist)=colnames(TRpw)
TRpw_dist
print(mn<-mean(TRpw_dist))
print(med<-median(TRpw_dist))
print(max<-max(TRpw_dist))
print(min<-min(TRpw_dist))

MLLtab_pt02=as.data.frame(cbind('sample'=row.names(genclone@tab),'MLL'=mll(genclone)))
MLLtab_pt02[MLLtab_pt02$sample %in% TRs,]
# this is calling only one pair of TRs (out of 6) a clonal group

# what if we go with 0.06209936 (max dist between TRs) as threshold (or a little above)
# note that pw distances arent exact thresholds because the calculated nodes are based on clusters!
# so set threshold above node value of TR clusters
mlg.filter(genclone, distance = prev_dist,algorithm = 'f') <- 0.08
genclone
# 99 MLGs

MLLtab_pt08=as.data.frame(cbind('sample'=row.names(genclone@tab),'MLL'=mll(genclone)))
MLLtab_pt08[MLLtab_pt08$sample %in% TRs,]

# now zoom in on histogram and look for elbow to find cutoff for pop gen
library(reshape2)
p=ggplot(melt(as.matrix(prev_dist)),aes(x=value))+
  geom_histogram(color='black',fill='lightgray')+
  theme_classic(base_size=14)
p+xlim(0,0.15)+geom_vline(xintercept=0.105)+
  geom_vline(xintercept=TRpw_dist,lty=3)+
  geom_vline(xintercept=farthest_thresh,color='red')
# will go with 0.105 to be conservative and exclude any highly related individuals
mlg.filter(genclone, distance = prev_dist,algorithm = 'f') <- 0.105
genclone
# 90 MLGs with 0.105 threshold (89 excluding CC7) -- this is what we are going to go with
MLLtab_pt105=as.data.frame(cbind('sample'=row.names(genclone@tab),'MLL'=mll(genclone)))
MLLtab_pt105[MLLtab_pt105$sample %in% TRs,]
# merge with coverage data, then can select one individual from clonal group with highest cov
cov=read.table('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2023/qranksAip.txt')
colnames(cov)=c('sample','sites5x')
MLLwCov=merge(MLLtab_pt105,cov,by='sample',all.x=T)
MLLpt02wCov=merge(MLLtab_pt02,cov,by='sample',all.x=T)
library(dplyr)
samples4popgen <- MLLwCov %>% 
  group_by(MLL) %>%
  filter(sites5x == max(sites5x)) %>%
  as.data.frame()
samples4popgen_pt02 <- MLLpt02wCov %>% 
  group_by(MLL) %>%
  filter(sites5x == max(sites5x)) %>%
  as.data.frame()
# export list of samples to use for pop gen
writeLines(samples4popgen$sample[!samples4popgen$sample=='CC7_host.bam'],
            '/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLLs89')
writeLines(samples4popgen_pt02$sample[!samples4popgen_pt02$sample %in% c(TR2s,'CC7_host.bam')],
           '/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLLs184')
# export sample to MLL tables
write.csv(MLLtab_pt02,'/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLL184tab_pt023.csv')
write.csv(MLLtab_pt08,'/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLL98tab_pt08.csv')
write.csv(MLLtab_pt105,'/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLL89tab_pt105.csv')

#### stats ####
# let's remove tech reps (+CC7) and run stats
mlg.filter(genclone, distance = prev_dist, algorithm = "f") <- farthest_thresh
genclone
keep=labels(prev_dist)[!labels(prev_dist) %in% c(TR2s,'CC7_host.bam')]
gensub=genclone[keep,]
gensub
sum=poppr(gensub,sample=99)
# cool 185 genets
#write.csv(sum,'/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/bams300_minusTRsCC7_poppr_stats.csv')
# are there clones across sites?
mlg.crosspop(gensub,strata=~POP)
mlg.table(gensub, strata=~POP)
# are there clones across subsites?
mlg.crosspop(gensub,strata=~SUB)
#mlg.table(gensub, strata=~SUB)
# are there clones across roots?
mlg.crosspop(gensub, strata=~RT)
#mlg.table(gensub, strata=~RT)

# psex
#(G <- rowSums(mlg.table(gensub, plot = FALSE) > 0))
ppsex <- psex(gensub, by_pop = FALSE, d = "rrmlg", mul = 1/2, method = "multiple")
quartz()
plot(ppsex, log = "y", col = ifelse(ppsex > 0.05, "red", "blue"))
dev.off()
# this isnt working

# then for pop gen
mlg.filter(genclone, distance = prev_dist, algorithm = "f") <- 0.105
genclone
keep=labels(prev_dist)[!labels(prev_dist) %in% c(TR2s,'CC7_host.bam')]
gensub=genclone[keep,]
gensub
sum_popogen=poppr(gensub,sample=99)
# 89 MLLs
write.csv(sum_popogen,'/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/popgen_pt105_bams300_minusTRsCC7_poppr_stats.csv')
tab=mlg.table(gensub,strata=~RT)
rowSums(tab>0)
# some NAs, need to exclude these in calculations
hist(rowSums(tab>0)[!grepl("NA",names(rowSums(tab>0)))])
mean(rowSums(tab>0)[!grepl("NA",names(rowSums(tab>0)))])
max(rowSums(tab>0)[!grepl("NA",names(rowSums(tab>0)))])
mean(rowSums(tab>0)[!grepl("NA|SI",names(rowSums(tab>0)))])
max(rowSums(tab>0)[!grepl("NA|SI",names(rowSums(tab>0)))])
genets_per_rt=as.data.frame(rowSums(tab>0)[!grepl("NA",names(rowSums(tab>0)))])
colnames(genets_per_rt)="n_genets"
write.csv(genets_per_rt,"/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/stats/n_MLLs10_per_rt.csv")

# now heterozygosity and inbreeding using hierfstat
library(hierfstat)
basicstat <- basic.stats(genclone, diploid = TRUE, digits = 2) 
basicstat$overall
Ho <- colMeans(basicstat$Ho,na.rm=T)
He <- colMeans(basicstat$Hs,na.rm=T)
Fis<-colMeans(basicstat$Fis,na.rm = T)
gen_stats <- as.data.frame(cbind(Ho,He,Fis))
gen_stats$Site=row.names(gen_stats)
write.csv(gen_stats,'/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/hierfstats_bams300.csv')

#### plots ####
# plot dendogram
hc=hclust(prev_dist,"complete") # "complete" method, is farthest neighbor, "ave" is UPGMA
par(mar = c(0, 2.2, 0, 0))
plot(hc,labels=F,ann=F)
abline(h=farthest_thresh,col="red")
#abline(h=TRpw_dist,col="black",lty=3)
#abline(h=0.08,col="green")

# plot histogram again with tech reps
quartz()
genclone_filtered <- filter_stats(genclone, distance = prev_dist, plot = TRUE)
abline(v=TRpw_dist,col="black",lty=3)
abline(v=farthest_thresh,col="red")
abline(v=0.105,col="black")
quartz.save('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/bams300_poppr_filterstats.pdf',type='pdf')
dev.off()

# plot histogram zoomed in
p=ggplot(melt(as.matrix(prev_dist)),aes(x=value))+
  geom_histogram(color='black',fill='lightgray')+
  theme_classic(base_size=14)
p+xlim(0,0.15)+geom_vline(xintercept=0.105)+
  geom_vline(xintercept=TRpw_dist,lty=3)+
  geom_vline(xintercept=farthest_thresh,color='red')+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank())+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

# color TRs
col_TRs=rep('black',length(labels(prev_dist)))
col_TRs[labels(prev_dist) %in% TRs] = 'red' 
library(sparcl)
quartz()
par(cex=1)
ColorDendrogram(hc, y = col_TRs, branchlength=0.05,ylab='prevosti')
abline(h=farthest_thresh,col="red")
#abline(h=TRpw_dist,col="black",lty=3)
abline(h=0.105,col="black")
#quartz.save('AipClean20_IBSpt17_wClones.pdf',type='pdf')
dev.off()

# color by site
site=sapply(strsplit(labels(prev_dist),'_'),FUN='[',1)
#siteCol=c("SI"="red","BP"="orange","LK"="yellow","TK"="green","MI"="blue","OK"="purple","CC7"="black")
siteCol2=readRDS('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/MS/final_figs/colPalSites.RDS')
siteCol=unlist(siteCol2$site)
TRlabs=labels(prev_dist)
TRlabs[!TRlabs %in% TRs]<-""
TRlabs[TRlabs %in% TRs]<-"*"
quartz()
par(cex=1,mar = c(1, 3, 1, 1))
ColorDendrogram(hc, y = siteCol[site], branchlength=0.055,ylab='prevosti distance',labels=TRlabs,xlab='')
abline(h=farthest_thresh,col="red")
abline(h=0.105,col="black")
#quartz.save('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/bams303_prev_dend_by_site.pdf',type='pdf')
dev.off()

# final -- more flexible plotting using dendextend
library(dendextend)
library(tidyverse)
library(stats)
dend <- as.dendrogram(hc)
pop.dend=sapply(strsplit(labels(dend),'_'),FUN='[',1)
#siteCol=c("SI"="red","BP"="orange","LK"="yellow","TK"="green","MI"="blue","OK"="purple","CC7"="black")
TRlabs=labels(dend)
TRlabs[!TRlabs %in% c(TRs)]<-""
TRlabs[TRlabs %in% TRs]<-"*"
dend2 <- assign_values_to_leaves_edgePar(dend=dend, value = siteCol[pop.dend] , edgePar = "col")
quartz()
dend2 %>% hang.dendrogram %>% set("branches_lwd", 1.3) %>% set("labels", TRlabs) %>% set("labels_cex",1.5) %>% plot(ylim=c(0,0.3))
abline(h=farthest_thresh,col="black",lty="dashed")
abline(h=0.105,col="black",lty="dashed")
quartz.save('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/MS/final_figs/bams300_prev_dend_by_site_20Jan2026.pdf',type='pdf')
dev.off()



