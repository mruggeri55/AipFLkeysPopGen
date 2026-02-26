setwd('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/pop_gen/')

#### PCoA of IBS matrix with no clones #####
# reading list of bam files = order of samples in IBS matrix
bams=read.table("bams90/bams90",header=F)[,1]
bams=sub("\\.bam","",bams,perl=T)

# read in IBS mat with no clones
ma_no_clones=as.matrix(read.table('bams90/bams90.ibsMat'))
dimnames(ma_no_clones)=list(bams,bams)

# create metadata table
site=sapply(strsplit(split='_',x=colnames(ma_no_clones)),FUN='[',1)
subsite=sapply(strsplit(split='_',x=colnames(ma_no_clones)),FUN='[',2)
root=sapply(strsplit(split='_',x=colnames(ma_no_clones)),FUN='[',3)
sample=sapply(strsplit(split='_',x=colnames(ma_no_clones)),FUN='[',4)
meta=as.data.frame(cbind(site,subsite,root,sample))
rownames(meta)=colnames(ma_no_clones)

# setting up colors for plotting
palette(rainbow(length(unique(site))))
colors=as.numeric(as.factor(site))
colpops=as.numeric(as.factor(sort(unique(site))))

# clustering / PCoA based on identity by state (IBS) based on single read resampling
# (for low and/or uneven coverage)
hc=hclust(as.dist(ma_no_clones),"ave")
quartz()
plot(hc,cex=0.5)  # this shows how similar unique MLGs are
dev.off()

col_clone=rep('black',nrow(ma_no_clones))
BP=bams[grep('BP_*',bams)]
SI=bams[grep('SI_*',bams)]
TK=bams[grep('TK_*',bams)]
LK=bams[grep('LK_*',bams)]
MI=bams[grep('MI_*',bams)]
OK=bams[grep('OK_*',bams)]


col_clone[bams %in% BP] = 'orange'
col_clone[bams %in% SI] = 'red'
col_clone[bams %in% TK] = 'green'
col_clone[bams %in% LK] = 'yellow'
col_clone[bams %in% MI] = 'blue'
col_clone[bams %in% OK] = 'purple'

quartz()
library(WGCNA)
library(sparcl)
ColorDendrogram(hclust(as.dist(ma_no_clones),"ave"), y = col_clone, branchlength=0.05)
dev.off()

#### now PCoA
library(vegan)
conds=data.frame(cbind(site))
pp0=capscale(ma_no_clones~1)
pp=capscale(ma_no_clones~site,conds)

# significance of by-site divergence
adonis(ma_no_clones~site,conds)

# eigenvectors
plot(pp0$CA$eig) 

axes2plot=c(1,2)  
#quartz()
library(adegenet) # for transp()
cmd=pp0  # change to cmd=pp to see constrained ordination (data projection to maximize by-site separation)
plot(cmd,choices=axes2plot,display="sites",type="n") # choices - axes to display
points(cmd,choices=axes2plot,pch=19,col=transp(colors,alpha=0.7))
#ordihull(cmd,choices= axes2plot,groups= conds$grp,draw="polygon",col=1+as.numeric(unique(as.factor(conds$grp))),label=T)
ordispider(cmd,choices= axes2plot,groups=conds$site,col="grey80")
#ordiellipse(cmd,choices= axes2plot,groups= conds$site,draw="polygon",col=colpops,label=T)

#quartz.save(file='MDSplot_bysite_noclones.pdf',type='pdf')
#dev.off()

#quartz()
# unscaled, to identify outliers
plot(cmd$CA$u[,axes2plot],pch=19,col=colors)
ordispider(cmd$CA$u[,axes2plot],groups=conds$site,col="grey80")
ordiellipse(cmd$CA$u[,axes2plot],groups= conds$site,draw="polygon",col=colpops,label=T)
#identify(cmd$CA$u[,axes2plot],labels=colnames(ma_no_clones),n=3,cex=0.7)

#quartz.save(file='MDSplot_bysite_noclones_unscaled.pdf',type='pdf')
#dev.off()

