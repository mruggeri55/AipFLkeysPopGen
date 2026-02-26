setwd("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/its2_type_profiles")
library(phyloseq)
#ps <- readRDS("ps.its2final.RDS")
ps<-readRDS("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/ps.its2.newMLL.RDS")

# Get Relative Abundance of syms
ps.rel <- transform_sample_counts(ps, function(x) x / sum(x))
ps2.rel=subset_samples(ps.rel,!name=='CC7') # remove CC7
tax=tax_table(ps2.rel)
colnames(tax)
tax2=tax[,c("Clade","Majority.ITS2.sequence","ITS2.type.profile")]
tax_table(ps2.rel)<-tax2
ps.rel=ps2.rel

type.rel <- tax_glom(ps.rel, taxrank = "Majority.ITS2.sequence" )
clad.rel <- tax_glom(ps.rel, taxrank = "Clade")

#### beta diversity ####
library(vegan)
bet.ps.prof <- betadisper(vegdist(ps.rel@otu_table),ps.rel@sam_data$site)
bet.ps.type <- betadisper(vegdist(type.rel@otu_table),type.rel@sam_data$site)
bet.ps.clade <- betadisper(vegdist(clad.rel@otu_table),clad.rel@sam_data$site)
# note default vegdist is bray-curtis

anova(bet.ps.prof)
anova(bet.ps.type)
anova(bet.ps.clade)

bet.ps.HSD.prof <- TukeyHSD(bet.ps.prof)
bet.ps.HSD.type <- TukeyHSD(bet.ps.type)
bet.ps.HSD.clad <- TukeyHSD(bet.ps.clade)

plot(bet.ps.HSD.prof,las=1)
plot(bet.ps.HSD.type,las=1)
plot(bet.ps.HSD.clad,las=1)
# type and clade level similar but profile different
# profile level mainly SI lower div than others
# type/clade level still signature of SI but also BP and LK lower diversity than others
# this makes sense because LK, BP, and SI dominated by A

#### adonis tests #####
# do sites have similar composition?
samdf=data.frame(ps.rel@sam_data)
adonis2(ps.rel@otu_table ~ site, data=samdf, permutations=999)
adonis2(type.rel@otu_table ~ site, data=samdf, permutations=999)
adonis2(clad.rel@otu_table ~ site, data=samdf, permutations=999)

library(funfuns)
pw_ad_prof=pairwise.adonis(ps.rel@otu_table, factors=samdf$site, permutations=999) 
pw_ad_type=pairwise.adonis(type.rel@otu_table, factors=samdf$site, permutations=999) 
pw_ad_clade=pairwise.adonis(clad.rel@otu_table, factors=samdf$site, permutations=999) 
# prof level -- every site significantly different
# type/clade level -- SI, BP, LK not different but all others are

# export tables
out='/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/stats/'
write.csv(pw_ad_prof,paste(out,'pw_ad_prof_bysite.csv',''))
write.csv(pw_ad_type,paste(out,'pw_ad_type_bysite.csv',''))
write.csv(pw_ad_clade,paste(out,'pw_ad_clade_bysite.csv',''))

write.csv(bet.ps.HSD.prof$group,paste(out,'bet_div_HSD_prof.csv',''))
write.csv(bet.ps.HSD.type$group,paste(out,'bet_div_HSD_type.csv',''))
write.csv(bet.ps.HSD.clad$group,paste(out,'bet_div_HSD_clade.csv',''))

# Fstat here represents ratio of between groups variance to within group variance
# can test whether this is correlated to over water distance

#### spatial autocorrelation ######
WaterDist=read.csv('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/gps/over_water_dist.csv')
WaterDist$PopPair=paste(WaterDist$Site1,WaterDist$Site2, sep = '_')
library(dplyr)
Dist = WaterDist %>% group_by(PopPair) %>% summarise(waterdist_km=mean(over_water_dist_km))
Dist$pop1=sapply(strsplit(Dist$PopPair,'_'),FUN='[',1)
Dist$pop2=sapply(strsplit(Dist$PopPair,'_'),FUN='[',2)
library(reshape2)
Distma=dcast(Dist[,-1],pop1~pop2,value.var = 'waterdist_km')[,-1]
rownames(Distma)=colnames(Distma)

pw_ad=pw_ad_prof ###### change to level of interest
pw_ad2=pw_ad %>% add_row(pairs = c("OK vs OK","MI vs MI","TK vs TK","LK vs LK","BP vs BP","SI vs SI"), F.Model = rep(0,times=6))
pw_ad$pop1=sapply(strsplit(pw_ad$pairs,' vs '),FUN='[',1)
pw_ad$pop2=sapply(strsplit(pw_ad$pairs,' vs '),FUN='[',2)
pw_ad2$pop1=sapply(strsplit(pw_ad2$pairs,' vs '),FUN='[',2)
pw_ad2$pop2=sapply(strsplit(pw_ad2$pairs,' vs '),FUN='[',1)
pw_ad=rbind(pw_ad,pw_ad2)
adma=dcast(pw_ad[,-1],pop1~pop2,value.var = 'F.Model')[,-1]
rownames(adma)=colnames(adma)

rownames(adma)==rownames(Distma)
library(ape)
mantel.test(adma,Distma)
# not significant spatial correlation at any level (prof,type, or clade)
# aka sites closer together do not have more similar communities
space=mantel(adma,Distma)

# what about temperature?
temp=read.csv("~/Desktop/Kenkel_lab/Aip_FLkeys/temp/buoy2019.csv")
sum=temp %>% group_by(site) %>% summarise(mean(Temp))
temp.dist=as.matrix(dist(sum,upper=T,diag=T))
colnames(temp.dist)<-sum$site
rownames(temp.dist)<-sum$site
# anova=aov(temp$Temp~temp$site)
# temp_pw=as.data.frame(TukeyHSD(anova)$`temp$site`)
# temp_pw2=temp_pw
# temp_pw$pop1=sapply(strsplit(rownames(temp_pw),'-'),FUN='[',1)
# temp_pw$pop2=sapply(strsplit(rownames(temp_pw),'-'),FUN='[',2)
# temp_pw2$pop1=sapply(strsplit(rownames(temp_pw2),'-'),FUN='[',2)
# temp_pw2$pop2=sapply(strsplit(rownames(temp_pw2),'-'),FUN='[',1)
# temp_pw=rbind(temp_pw,temp_pw2)
# temp_pw$diff=abs(temp_pw$diff)
#library(reshape2)
#tempma=dcast(temp_pw,pop1~pop2,value.var = 'diff')[,-1]
#rownames(tempma)=colnames(tempma)
#diag(tempma)<-0
#temp.dist=tempma

# match up buoy and site order
tempma=temp.dist[match(c("venice","virginia","BN","PK","vaca","key west"),rownames(temp.dist)),match(c("venice","virginia","BN","PK","vaca","key west"),colnames(temp.dist))]
adma=adma[match(c("OK","MI","TK","LK","BP","SI"),rownames(adma)),match(c("OK","MI","TK","LK","BP","SI"),colnames(adma))]
#temp.man=mantel.test(adma,tempma,nperm=9999,graph=T)
env=mantel(adma,tempma)
# clade significant!! r=0.6467,p=0.029
# prof ns

# adding in genetic correlation using Fst
Fst=read.csv("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/pop_gen/MLLs89/Fst_mat.csv",row.names = 1)
rownames(Fst)==rownames(adma)
host=mantel(adma,Fst)
# for clade, marginally significant (r=0.4221,p=0.06)

space
env
host
# together mantel tests at site level support that 
# higher order symbiont community structure (clade level)
# is defined by temperature rather than host genetic structure
# type also follows this pattern
# but prof ns across the board (note could maybe subset within each clade)
library(corrplot)
corrplot(as.matrix(Fst),type='lower', is.corr = F,diag = F,addCoef.col = 'black',
         tl.col = 'black', tl.cex=1.5,
         col.lim=c(0,0.1),title="host Fst")
corrplot(as.matrix(adma),type='lower', is.corr = F,diag = F,addCoef.col = 'black',
         tl.col = 'black', tl.cex=1.5,title="community dissimilarity")
corrplot(as.matrix(tempma),type='lower', is.corr = F,diag = F,addCoef.col = 'black',
         tl.col = 'black', tl.cex=1.5,title="temperature")
Distma<-Distma[match(c('OK','MI','TK','LK','BP','SI'),rownames(Distma)),match(c('OK','MI','TK','LK','BP','SI'),colnames(Distma))]
corrplot(as.matrix(Distma),type='lower', is.corr = F,diag = F,addCoef.col = 'black',
         tl.col = 'black', tl.cex=1.5,title="spatial distance")


#### RDA #####
comm<-ps.rel@otu_table ##### change to level of interest
sum$site<-c('TK','LK','SI','BP','OK','MI')
colnames(sum)[2]="temp"
env<-merge(samdf,sum,by="site",all.x=T,sort=F)[c("sampleID","site","temp")]
gps<-read.csv('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/gps/Site_GPS.csv')
colnames(gps)<-c("site","lat","long")
env<-merge(env,gps,by='site',all.x=T,sort=F)
rownames(env)<-env$sampleID
env$sampleID<-NULL
env<-env[match(rownames(comm),rownames(env)),]

res<-capscale(comm~.,data=env[,2:3,drop=F],distance = "bray")
summary(res)
anova(res)
plot(res)


##### spatial correlation of each profile using Morans I #####
Distma<-read.csv("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/spatial_distance_matrix_allsamps.csv",row.names = 1)
sym.abund<-otu_table(ps.rel) ##### change to level of interest
rownames(sym.abund)<-samdf$name
sym.abund.noNA<-sym.abund[match(rownames(Distma),rownames(sym.abund)),]
rownames(Distma)==rownames(sym.abund.noNA)
zoox_ma<-as.data.frame(sym.abund.noNA)

library(vegan)
library(ape)
Mtest <-data.frame()
dist.mat <- as.matrix(1/Distma) # inverse because closer observations should have higher weights
diag(dist.mat) <- 0 
for (col in colnames(zoox_ma[,colSums(zoox_ma>0)/nrow(zoox_ma)>0.05])) {
  x <- zoox_ma[,col]
  I <- Moran.I(x, dist.mat)
  Mtest[col,1]<-print(I$p.value)
  Mtest[col,2]<-print(I$observed)
  Mtest[col,3]<-print(I$expected)
  Mtest[col,4]<-print(I$sd)
}
colnames(Mtest)=c("p.val","I.obs","I.exp","I.sd")
Mtest$p.adj <- Mtest$p.val %>% p.adjust(method = "BH")
length(Mtest$I.exp[Mtest$p.adj<0.05])/length(Mtest$I.exp)*100
Mtest.spatial<-Mtest
# pretty strong spatial autocorrelation in syms
# now check if this is still there after removing putative clones
MLLs89=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLLs89")
MLLs89=gsub("_host.bam","",MLLs89$V1)
MLLs89[!MLLs89 %in% rownames(Distma)] # 13 samples not in distma
MLLs89sub=MLLs89[MLLs89 %in% rownames(Distma)]
Distma=as.matrix(Distma)
dist=Distma[match(MLLs89sub,rownames(Distma)),match(MLLs89sub,colnames(Distma))]
zoox_ma=zoox_ma[match(MLLs89sub,rownames(zoox_ma)),]

Mtest <-data.frame()
dist.mat <- as.matrix(1/dist) # inverse because closer observations should have higher weights
diag(dist.mat) <- 0 
for (col in colnames(zoox_ma[,colSums(zoox_ma>0)/nrow(zoox_ma)>0.05])) {
  x <- zoox_ma[,col]
  I <- Moran.I(x, dist.mat)
  Mtest[col,1]<-print(I$p.value)
  Mtest[col,2]<-print(I$observed)
  Mtest[col,3]<-print(I$expected)
  Mtest[col,4]<-print(I$sd)
}
colnames(Mtest)=c("p.val","I.obs","I.exp","I.sd")
Mtest$p.adj <- Mtest$p.val %>% p.adjust(method = "BH")
length(Mtest$I.exp[Mtest$p.adj<0.05])/length(Mtest$I.exp)*100
M.test.space.89<-Mtest

# export results
write.csv(Mtest.spatial,paste(out,'MoransI_allinds.csv'))
write.csv(M.test.space.89,paste(out,'MoransI_noclones.csv'))

#### note the above part randomly subsetted 89 individuals (aka not representing full community comp within clonal groups)
