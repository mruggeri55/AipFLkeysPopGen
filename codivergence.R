setwd('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025')

# read in sym rel abund and meta
zooxCounts = read.csv('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/its2_type_profiles/relAbundfromPhyloSeq.csv', row.names = 1)
zooxCounts$sample=rownames(zooxCounts)
meta_MLL=read.csv("syms_among_clones/meta_MLLs184.csv")
meta_MLL=meta_MLL[!grepl("TR2",meta_MLL$name),]
meta_MLL$MLL=as.factor(meta_MLL$MLL)
zooxCounts_master=merge(zooxCounts,meta_MLL,by='sample',all.x=T)
zooxCounts_master$sample[duplicated(zooxCounts_master$sample)]
rownames(zooxCounts_master)=zooxCounts_master$name

# read in prevosti distance matrix
prev_dist=readRDS("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/calling_clones/bams300_prev_dist.RDS")
samp_order=gsub('_host.bam',"",labels(prev_dist))

# reorder zoox df to match prev order
# also subset to only samples in both datasets
zoox_match=zooxCounts_master[rownames(zooxCounts_master) %in% samp_order,]
samp_order=samp_order[samp_order %in% rownames(zoox_match)]
zoox_order=zoox_match[order(rownames(zoox_match),samp_order),]

# now make into matrixes
prev_ma=as.matrix(prev_dist)
rownames(prev_ma)=gsub("_host.bam","",rownames(prev_ma))
colnames(prev_ma)=gsub("_host.bam","",colnames(prev_ma))
prev_sub=prev_ma[rownames(prev_ma) %in% samp_order,
                 colnames(prev_ma) %in% samp_order]
nrow(prev_sub)
prev_ord=prev_sub[order(rownames(prev_sub),samp_order),
                 order(colnames(prev_sub), samp_order)]

colnames(zoox_order)
#zoox_ma=as.matrix(zoox_order[,2:36])
rownames(zoox_order)==rownames(prev_ord)
colnames(zoox_order)=gsub('\\.','-',colnames(zoox_order))
colnames(zoox_order)
colnames(zoox_order)[c(3,9,11,29,30,31,32,35)]=c("A4/A4m","A1/A4-A1ac-A1bn","A4/A1","C3/C1","C15/C116aa", "C3/C15" ,"C15/C15ch","D1/D2")
colnames(zoox_order)
zoox_ma=as.matrix(zoox_order[,2:36])

write.csv(zoox_ma,"/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/zoox_ma_Jul2025.csv")

##### plotting #####
library(pheatmap)
sites=as.data.frame(sapply(strsplit(x=rownames(prev_ord),split='_'),FUN='[',1))
rownames(sites)=rownames(prev_ord)
colnames(sites)='site'
hc=hclust(as.dist(prev_ord),"ave")

library(RColorBrewer)
cols <- c("white",colorRampPalette(brewer.pal(9, "Blues"))(25))
ann_colors=list(site=c(SI='#FF6663',BP='#FEB144',LK='#FDFD97',TK='#9EE09E',MI='#9EC1CF',OK='#CC99C9',CC7='black'))

pheatmap(t(zoox_ma),cluster_cols = hc, cluster_rows = F, show_colnames = F,
         col=cols,annotation_col = sites, annotation_colors=ann_colors)

# now let's try to add tree for sym profs, subset each genus
SymProfDist=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_profile_distances/A/20220808T040613_unifrac_profile_distances_A_sqrt.dist",row.names = 1)[,-1]
BrevProfDist=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_profile_distances/B/20220808T040613_unifrac_profile_distances_B_sqrt.dist",row.names = 1)[,-1]
CladProfDist=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_profile_distances/C/20220808T040613_unifrac_profile_distances_C_sqrt.dist",row.names = 1)[,-1]
DurProfDist=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_profile_distances/D/20220808T040613_unifrac_profile_distances_D_sqrt.dist",row.names = 1)[,-1]

colnames(SymProfDist)=rownames(SymProfDist)
colnames(BrevProfDist)=rownames(BrevProfDist)
colnames(CladProfDist)=rownames(CladProfDist)
colnames(DurProfDist)=rownames(DurProfDist)

symHC=hclust(as.dist(SymProfDist),"ave")
bHC=hclust(as.dist(BrevProfDist),"ave")
cHC=hclust(as.dist(CladProfDist),"ave")
dHC=hclust(as.dist(DurProfDist),"ave")

symA=as.matrix(zoox_ma[,grep("A",colnames(zoox_ma))])
symB=as.matrix(zoox_ma[,grep("B",colnames(zoox_ma))])
symC=as.matrix(zoox_ma[,grep("C",colnames(zoox_ma))])
symD=as.matrix(zoox_ma[,grep("D",colnames(zoox_ma))])

# re-order profiles in sym matrix to match cluster order
# do not change row order since this matches with host relatedness
symAord=symA[,match(colnames(SymProfDist),colnames(symA))]
symBord=symB[,match(colnames(BrevProfDist),colnames(symB))]
symCord=symC[,match(colnames(CladProfDist),colnames(symC))]
symDord=symD[,match(colnames(DurProfDist),colnames(symD))]

# plot heatmaps
pA=pheatmap(t(symAord),cluster_cols = hc, cluster_rows = symHC, 
            show_colnames = F,color=cols,annotation_col = sites,
            annotation_colors=ann_colors,breaks=seq(from=0,to=1,by=0.04),
            treeheight_col = 70,legend=F,annotation_legend = F,annotation_names_col=F,
            labels_row = rep('',length(colnames(symAord))),cutree_cols = 90)
pB=pheatmap(t(symBord),cluster_cols = hc, cluster_rows = bHC, 
            show_colnames = F,color=cols,breaks=seq(from=0,to=1,by=0.04),
            treeheight_col = 0,legend=F,labels_row = rep('',length(colnames(symBord))),cutree_cols = 90)
pC=pheatmap(t(symCord),cluster_cols = hc, cluster_rows = cHC, 
            show_colnames = F,color=cols,breaks=seq(from=0,to=1,by=0.04),
            treeheight_col = 0,legend=F,labels_row = rep('',length(colnames(symCord))),cutree_cols = 90)
pD=pheatmap(t(symDord),cluster_cols = hc, cluster_rows = dHC, 
            show_colnames = F,color=cols,breaks=seq(from=0,to=1,by=0.04),
            treeheight_col = 0,legend=F,labels_row = rep('',length(colnames(symDord))),cutree_cols = 90)

library(cowplot)
library(ggpubr)
library(gridExtra)
#quartz()
plot_grid(pA[[4]],pB[[4]],pC[[4]],pD[[4]],ncol=1,align='v',axis='lr',rel_heights = c(25,10,9,3))
#dev.off()

##### check whether symbiont communities are correlated to host genetics ####
# subsetting one sample per clonal group using smallest threshold (2.3%)
# doing this randomly, will have to think about combining communities within clonal groups
MLLs184sym=zoox_order[unique(zoox_order$MLL),2:36]
MLLs184host=prev_ord[rownames(prev_ord) %in% rownames(MLLs184sym),colnames(prev_ord) %in% rownames(MLLs184sym)]
# subsetting only unique genets using greatest threshold (10%)
MLLs89=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLLs89")
MLLs89sym=zoox_order[zoox_order$name %in% gsub("_host.bam","",MLLs89$V1),2:36]
MLLs89host=prev_ord[rownames(prev_ord) %in% rownames(MLLs89sym),colnames(prev_ord) %in% rownames(MLLs89sym)]

##### try between sample distances from symportal ####
sampA=as.matrix(read.delim("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_sample_distances/A/20220808T040613_unifrac_sample_distances_A_sqrt.dist",header=F,row.names = 1)[,-1])
sampB=as.matrix(read.delim("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_sample_distances/B/20220808T040613_unifrac_sample_distances_B_sqrt.dist",header=F,row.names = 1)[,-1])
sampC=as.matrix(read.delim("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_sample_distances/C/20220808T040613_unifrac_sample_distances_C_sqrt.dist",header=F,row.names = 1)[,-1])

# change the following to run each sym
clade="Symbiodinium"
sym.dist=sampA

# rownames are sample names
# second column is symportal ID so do not read in
colnames(sym.dist)=rownames(sym.dist)
sym.distsub=sym.dist[rownames(sym.dist) %in% zoox_order$sample,colnames(sym.dist) %in% zoox_order$sample]
# note some samples will be missing if they do not have any A syms!
samp_list=as.data.frame(cbind("name"=zoox_order$name[zoox_order$sample %in% rownames(sym.distsub)],"sample"=zoox_order$sample[zoox_order$sample %in% rownames(sym.distsub)]))
prevAsub=as.matrix(prev_ord[rownames(prev_ord) %in% samp_list$name,colnames(prev_ord) %in% samp_list$name])
prevAord=as.matrix(prevAsub[match(samp_list$name,rownames(prevAsub)),match(samp_list$name,colnames(prevAsub))])
sym.distord=sym.distsub[match(samp_list$sample,rownames(sym.distsub)),match(samp_list$sample,colnames(sym.distsub))]

library(vegan)
library(ape)
mantel.test(sym.distord,prevAord)
correlog <- mantel.correlog(sym.distord,prevAord, nperm = 999) # nperm for permutations
plot(correlog)

samp_list_MLL89=samp_list[samp_list$name %in% gsub("_host.bam","",MLLs89$V1),]
mantel.test(sym.distord[rownames(sym.distord) %in% samp_list_MLL89$sample,colnames(sym.distord) %in% samp_list_MLL89$sample],
            prevAord[rownames(prevAord) %in% samp_list_MLL89$name,colnames(prevAord) %in% samp_list_MLL89$name])
# getting significant correlation for both of these
# but looking at plots below may be heteroskedastic
library(reshape2)
long=cbind("sym.dist"=melt(sym.distord)[,3],melt(prevAord, value.name = "host.dist"))
library(ggplot2)
ggplot(long[!long$host.dist==0,],aes(x=sym.dist,y=host.dist))+
  geom_point(size=0.2)+geom_smooth(method="lm")+
  ggtitle(clade)
# only unique genotypes
ggplot(long[long$Var1 %in% samp_list_MLL89$name & long$Var2 %in% samp_list_MLL89$name & !long$host.dist==0,],
       aes(x=sym.dist,y=host.dist))+geom_point(size=0.2)+geom_smooth(method="lm")
# also removing SI
ggplot(long[long$Var1 %in% samp_list_MLL89$name & long$Var2 %in% samp_list_MLL89$name & !long$host.dist==0 & !grepl("SI",long$Var1) & !grepl("SI",long$Var2),],
       aes(x=sym.dist,y=host.dist))+
  geom_point(size=0.2)+geom_smooth(method="lm")
samp_list_MLL89noSI=samp_list[samp_list$name %in% gsub("_host.bam","",MLLs89$V1) & !grepl("SI",samp_list$name),]
mantel.test(sym.distord[rownames(sym.distord) %in% samp_list_MLL89noSI$sample,colnames(sym.distord) %in% samp_list_MLL89noSI$sample],
            prevAord[rownames(prevAord) %in% samp_list_MLL89noSI$name,colnames(prevAord) %in% samp_list_MLL89noSI$name])
# not sig without clones or SIs

# checking heteroskadascity
lm=lm(host.dist~sym.dist,data=long[long$Var1 %in% samp_list_MLL89$name & long$Var2 %in% samp_list_MLL89$name & !long$host.dist==0,])
qqPlot(residuals(lm))
plot(lm,which=1)
# very bad with clones but not terrible without clones
# maybe can normalize prev dist? squaring seems to help, note sym dist already on sqrt scale
lm=lm(host.dist^2~sym.dist,data=long[!long$host.dist==0,])
qqPlot(residuals(lm))
plot(lm,which=1)


#### Moran test ####
# following https://github.com/baumlab/Starko_et_al_Porites_KI/blob/main/Analyses/CoevolutionScript.R
Mtest <-data.frame()
dist.mat <- 1/prev_ord
diag(dist.mat) <- 0
rownames(zoox_ma)==rownames(prev_ord)
# note "C15.C116aa" and "B1.B1g" all zeros (must be from subsetting)
# exclude these columns for loop to work
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
M.test.codiv.all<-Mtest
# all significant associations (note including clones) 
# now without clones
Mtest2=data.frame()
dist.mat2=1/MLLs89host
diag(dist.mat2) <- 0
for (col in colnames(MLLs89sym[,colSums(MLLs89sym>0)/nrow(MLLs89sym)>0.05])) {
  x <- MLLs89sym[,col]
  I <- Moran.I(x, dist.mat2)
  Mtest2[col,1]<-print(I$p.value)
  Mtest2[col,2]<-print(I$observed)
  Mtest2[col,3]<-print(I$expected)
  Mtest2[col,4]<-print(I$sd)
}
colnames(Mtest2)=c("p.val","I.obs","I.exp","I.sd")
Mtest2$p.adj <- Mtest2$p.val %>% p.adjust(method = "BH")
length(Mtest2$I.exp[Mtest2$p.adj<0.05])/length(Mtest2$I.exp)*100
# only A profs sig after removing clones
M.test.codiv.89<-Mtest2

#### try Paco ####
library(paco)
prefix="A_min5ind_noSIs"
sub=zoox_ma[,grep("A",colnames(zoox_ma))] # for all
#sub=zoox_ma[rownames(zoox_ma) %in% gsub("_host.bam","",MLLs89$V1), 
#            grep("A",colnames(zoox_ma))] # for no clones

sub2=sub[rowSums(sub>0)>=1,colSums(sub>0)>=5] # subset only variants present in 5 or more samples
sub2=sub2[!grepl("SI",rownames(sub2)),]
#sub=sub[,match(colnames(SymProfDist),colnames(sub))]
#colnames(sub)==colnames(SymProfDist)
sub=sub2
ProfDist=as.matrix(SymProfDist[match(colnames(sub2),rownames(SymProfDist)),match(colnames(sub2),colnames(SymProfDist))]) # change to match genus
#HostDist=prev_ord # for all
#HostDist=prev_ord[rownames(prev_ord) %in% gsub("_host.bam","",MLLs89$V1),colnames(prev_ord) %in% gsub("_host.bam","",MLLs89$V1)] # for no clones
HostDist=prev_ord[match(rownames(sub),rownames(prev_ord)),match(rownames(sub),colnames(prev_ord))]
rownames(HostDist)==rownames(sub)

paco_data <- prepare_paco_data(H = HostDist, P = ProfDist, sub)
paco_data <- add_pcoord(paco_data, correction = 'cailliez')
z_final <- PACo(paco_data, nperm = 1000, seed = 999, symmetric = FALSE, shuffled = TRUE)
# https://cran.r-project.org/web/packages/paco/paco.pdf
z_final$gof

hist(z_final$shuffled, xlim=c(z_final$gof$ss-1, max(z_final$shuffled)+1), las = 1, col = "lightblue", xlab = "Sum of squared residuals", main = '')
abline(v = z_final$gof$ss, col = "blue", lwd = 3)

# following supp file 1 of below
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061048#s5
HP.proc <- procrustes(z_final$H_PCo,z_final$P_PCo)
# HostX <- HP.proc$X
# ParY <- HP.proc$Yrot
# plot(HostX, asp=1, pch=16,col="blue")
# arrows(ParY[,1], ParY[,2], HostX[,1], HostX[,2], length=0.12, angle=15,
#        xpd=FALSE) 
# points(ParY, pch=2,col="red")
# identify(ParY[,1], ParY[,2], rownames(ParY), offset=0.3, xpd=FALSE, cex=0.8)
# # click on sym points to display
# identify(HostX[,1], HostX[,2], rownames(HostX),offset=0.3, xpd=TRUE, cex= 0.8)
# then do the same with host points

df=as.data.frame(HP.proc$X)
df$site=as.factor(sapply(strsplit(rownames(df),'_'),FUN='[',1))
ggplot(df,aes(x=Axis.1,y=Axis.2,color=site))+
  geom_point()+
  scale_color_manual(values=ann_colors$site)+
  theme_bw()

symdf=as.data.frame(HP.proc$Yro)
symdf$site=gsub("[.]\\d+$", "",rownames(symdf))
colnames(symdf)=colnames(df)
symdf$shape=rep("sym",length(symdf$Axis.1))
# site is really sym but need same naming so can color properly
ggplot(symdf,aes(x=Axis.1,y=Axis.2,color=site))+  
  geom_point()

df$shape=rep("host",length(df$Axis.1))
all=rbind(df,symdf)
ggplot(all,aes(x=Axis.1,y=Axis.2,color=site,shape=shape))+
  geom_point(size=2)+
  scale_color_manual(values=ann_colors$site)+
  theme_bw()
# don't know how to get arrows on this
# maybe use different columns?
symdf$sym=gsub("[.]\\d+$", "",rownames(symdf))
test=cbind(df[,c("Axis.1","Axis.2","site")],symdf[,c("Axis.1","Axis.2","sym")])
colnames(test)[4:5]=c("sym.Axis.1","sym.Axis.2")
#A4 = colorRampPalette(c("navy", "turquoise3"))
#A = as.list(A4(10))
#names(A)=unique(test$sym)
#ann_syms=list("A4"="#000080", "A4.A4m"="#001588", "A4.A4cs" ="#002B91","A4.A4cu.A4ab.A4t.A4cv"= "#004199","A4.A4q.A4cw"  = "#0057A2","A4.A4m.A4cs.A4ct" = "#006DAA","A4.A4ct.A4cs.A4q" = "#0083B3","A4.A4q.A1il"="#0099BB" ,"A4.A4ct.A4cs.A1il" = "#00AFC4" ,"A4.A1"= "#00C5CD")
ann_syms=list("A4"="#000080", "A4.A4m"="#001588", "A4.A4cs" ="#002B91","A4.A4q.A4cw"  = "#0057A2",
              "A4.A4m.A4cs.A4ct" = "#006DAA","A4.A4ct.A4cs.A4q" = "#0083B3","A4.A4q.A1il"="#0099BB" ,"A4.A4ct.A4cs.A1il" = "#00AFC4","A1.A4.A1ac.A1bn"="#00C5CD")
#ann_syms=list("A4"="#000080","A4.A4m"  = "#006DAA","A4.A4ct.A4cs.A4q" = "#0099BB","A4.A4q.A1il"="#00C5CD")

ggplot(test)+
  geom_segment(aes(x=Axis.1,xend=sym.Axis.1,y=Axis.2,yend=sym.Axis.2,color=sym))+
  geom_point(aes(x=Axis.1,y=Axis.2,color=site),size=2)+
  scale_color_manual(values=c(ann_syms,ann_colors$site))+
  geom_point(aes(x=sym.Axis.1,sym.Axis.2,color=sym,size=2),shape=17)+
  theme_bw()
  #facet_wrap(~sym)


saveRDS(z_final,paste("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/",
"paco_",prefix,".RDS",sep=""))

# A: co-phylogeny significant with or without clones (MLLs89)
# B: co-phylogeny sig with clones, ns without clones (MLLs89, note even more reduced sample size bc only counts those with B present)
# C: co-phylogeny ns without clones (MLLs89, note even more reduced sample size bc only counts those with C present)

##### trying tanglegram
#z_final=readRDS("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/paco_A_MLLs89.RDS")
library(ape)
library(dendextend)
library(phytools)
# Convert trees to dendrograms
host_dend <- as.dendrogram(hclust(as.dist(z_final$H)))
para_dend <- as.dendrogram(hclust(as.dist(z_final$P)))
assoc <- z_final$HP
assoc_mat <- assoc[labels(host_dend), labels(para_dend)]  
tanglegram(host_dend, para_dend,
           lab.cex = 0.8,
           lwd = 1,
           color_lines = assoc_mat,  # custom line coloring using matrix
           main = "PACo Host–SymA Tanglegram",
           common_subtrees_color_lines = FALSE)
host.tree<-as.phylo(host_dend)
sym.tree<-as.phylo(para_dend)
cophyloplot(host.tree, sym.tree, assoc = assoc,show.tip.label = F)
library(Rtapas)
LFc <- max_cong(assoc[host.tree$tip.label,sym.tree$tip.label], host.tree, sym.tree, n=5, N=1e4, method = "paco", symmetric = TRUE,
                ei.correct = "sqrt.D", percentile = 0.01, res.fq = TRUE)
x11()
tangle_gram(host.tree,sym.tree,HS=assoc[host.tree$tip.label,sym.tree$tip.label],LFc,
            method="paco",colscale = "sequential", colgrad = c("darkblue","darkblue"),
            node.tag = F,link.lwd=1.5)
dev.print(pdf,"/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/A_min5inds_Tanglegram_Jul2025.pdf") 
dev.off()

#### paco final ####
##### adding in MLL info to ps, then merging MLLs, then cophylogeny #####
library(phyloseq)
# ps <- readRDS("ps.its2final.RDS")
# samptoMLLs89=read.csv("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLL89tab_pt105.csv",row.names=1)
# colnames(samptoMLLs89)=c("name","MLLs89")
# samptoMLLs89$name=gsub("_host.bam","",samptoMLLs89$name)
# samptoMLLs184=read.csv("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLL184tab_pt023.csv",row.names=1)
# colnames(samptoMLLs184)=c("name","MLLs184")
# samptoMLLs184$name=gsub("_host.bam","",samptoMLLs184$name)
# samdf=data.frame(ps@sam_data) # extract metadata info
# samdf=samdf[,!colnames(samdf) %in% c("genotype","site_genet","sample")] # remove old genotype designations and redundant sample col
# samdfnew=merge(samdf,samptoMLLs184,by="name",all.x=T)
# samdfnew2=merge(samdfnew,samptoMLLs89,by="name",all.x=T)
# rownames(samdfnew2)<-samdfnew2$sampleID
# samdfnew2$site_MLLs89=paste(samdfnew2$site,samdfnew2$MLLs89,sep='-')
# samdfnew2$site_MLLs184=paste(samdfnew2$site,samdfnew2$MLLs184,sep='-')
# sample_data(ps) <- sample_data(samdfnew2) #change sample data to new df
# rownames(ps@otu_table)==rownames(ps@sam_data) # check if they are in the same order
# rownames(ps@otu_table)==ps@sam_data$sampleID
# tax=tax_table(ps)
# colnames(tax)
# tax2=tax[,c("Clade","Majority.ITS2.sequence","ITS2.type.profile")]
# tax_table(ps)<-tax2
# # save new phyloseq object with MLL designations
#saveRDS(ps,"/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/ps.its2.newMLL.RDS")

ps<-readRDS("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/ps.its2.newMLL.RDS")
MLLs89=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLLs89")
prev_dist=readRDS("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/calling_clones/bams300_prev_dist.RDS")
prev_ma=as.matrix(prev_dist)

SymProfDist=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_profile_distances/A/20220808T040613_unifrac_profile_distances_A_sqrt.dist",row.names = 1)[,-1]
BrevProfDist=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_profile_distances/B/20220808T040613_unifrac_profile_distances_B_sqrt.dist",row.names = 1)[,-1]
CladProfDist=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_profile_distances/C/20220808T040613_unifrac_profile_distances_C_sqrt.dist",row.names = 1)[,-1]
DurProfDist=read.table("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/between_profile_distances/D/20220808T040613_unifrac_profile_distances_D_sqrt.dist",row.names = 1)[,-1]
colnames(SymProfDist)=rownames(SymProfDist)
colnames(BrevProfDist)=rownames(BrevProfDist)
colnames(CladProfDist)=rownames(CladProfDist)
colnames(DurProfDist)=rownames(DurProfDist)

#ps.MLLs184 <- merge_samples(ps, "site_MLLs184") # collapse by MLLs
ps.MLLs89 <- merge_samples(ps, "site_MLLs89") # collapse by MLLs

# convert otu_table from abundance to presence absence
ps.MLLs89.bin <- transform_sample_counts(ps.MLLs89, function(x) {
  ifelse(x > 0, 1, 0)
})
ps.MLLs89.bin <- subset_samples(ps.MLLs89.bin, !is.na(MLLs89)) # remove samples without host genet info
# then need to convert sample names in MLLs184 host matrix to MLL to match
head(MLLs89$V1)
prev.89=prev_ma[rownames(prev_ma) %in% MLLs89$V1,
                 colnames(prev_ma) %in% MLLs89$V1]

samptoMLLs89=read.csv("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLL89tab_pt105.csv",row.names=1)
sub=samptoMLLs89[samptoMLLs89$sample %in% colnames(prev.89),]
sub$site_MLL=paste(sapply(strsplit(sub$sample,"_"),FUN='[',1),
                             sub$MLL,sep='-')
colnames(prev.89)==sub$sample
rownames(ps.MLLs89.bin@otu_table)==sub$site_MLL
# need to reorder tables to match otu table order
# rename host table with MLL ID
colnames(prev.89)<-sub$site_MLL
rownames(prev.89)<-sub$site_MLL
# now match host table cols and rows to otu table order
prev.89.2=prev.89[match(rownames(ps.MLLs89.bin@otu_table),rownames(prev.89)),
                       match(rownames(ps.MLLs89.bin@otu_table),colnames(prev.89))]
rownames(ps.MLLs89.bin@otu_table)==rownames(prev.89.2)
badMLLs=rownames(ps.MLLs89.bin@otu_table)[!rownames(ps.MLLs89.bin@otu_table) %in% rownames(prev.89)]
# this does not make sense, some genets in otu table not in prevosti table...
# that is not possible because they would not have genet designations if host wasn't sequenced
# figured this out, basically this an issue with the sample with highest host cov which got chosen for pop gen did not have its symbiont community sequenced!
# "BP-18"  "BP-242" "CC7-83" "LK-16"  "TK-239" -- could probably correct this for the BPs because multiple clones but will move forward for now
ps.MLLs89.bin<-prune_samples(!sample_names(ps.MLLs89.bin) %in% badMLLs,ps.MLLs89.bin)
prev.89.2=prev.89[match(rownames(ps.MLLs89.bin@otu_table),rownames(prev.89)),
                  match(rownames(ps.MLLs89.bin@otu_table),colnames(prev.89))]
rownames(ps.MLLs89.bin@otu_table)==rownames(prev.89.2)
# FINALLY THESE MATCH -- now can do codivergence
library(paco)
prefix="B_MLLs89comb"
ps.sub=subset_taxa(ps.MLLs89.bin,Clade=="B") #change to clade of interest
sub=as.matrix(ps.sub@otu_table)
#sub2=sub[,colSums(sub>0)>=5] # subset only variants present in 5 or more samples
#sub=sub2
ProfDist=BrevProfDist # change to match genus
ProfDist=as.matrix(ProfDist[match(colnames(sub),rownames(ProfDist)),
                               match(colnames(sub),colnames(ProfDist))]) 
HostDist=prev.89.2
rownames(HostDist)==rownames(sub)
# cannot do B with frequency of 5 inds because only leaves 2 profs

paco_data <- prepare_paco_data(H = HostDist, P = ProfDist, sub)
paco_data <- add_pcoord(paco_data, correction = 'cailliez')
z_final <- PACo(paco_data, nperm = 1000, seed = 999, symmetric = FALSE, shuffled = TRUE,method="r0")
# https://cran.r-project.org/web/packages/paco/paco.pdf
# note there are default methods! use method="" to set
# default method (r0) assumes syms track host
# c0 assumes host tracks syms
# swap/bakctrack does not specify
# so far for As -- r0, swap and backtrack significant but c0 not! 
saveRDS(z_final,paste("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/",
                      "paco_",prefix,".RDS",sep=""))

z_final$gof

hist(z_final$shuffled, xlim=c(z_final$gof$ss-1, max(z_final$shuffled)+1), las = 1, col = "lightblue", xlab = "Sum of squared residuals", main = '')
abline(v = z_final$gof$ss, col = "blue", lwd = 3)

z_final<-paco_links(z_final)


HP.proc <- procrustes(z_final$H_PCo,z_final$P_PCo)
df=as.data.frame(HP.proc$X)
df$site=as.factor(sapply(strsplit(rownames(df),'[.]'),FUN='[',1))
symdf=as.data.frame(HP.proc$Yro)
symdf$site=gsub("[.]\\d+$", "",rownames(symdf))
colnames(symdf)=colnames(df)
symdf$shape=rep("sym",length(symdf$Axis.1))
df$shape=rep("host",length(df$Axis.1))
all=rbind(df,symdf)
symdf$sym=gsub("[.]\\d+$", "",rownames(symdf))
test=cbind(df[,c("Axis.1","Axis.2","site")],symdf[,c("Axis.1","Axis.2","sym")])
colnames(test)[4:5]=c("sym.Axis.1","sym.Axis.2")
#ann_syms=list("A4"="#000080", "A4.A4m"="#001588", "A4.A4cs" ="#002B91","A4.A4q.A4cw"  = "#0057A2",
#              "A4.A4m.A4cs.A4ct" = "#006DAA","A4.A4ct.A4cs.A4q" = "#0083B3","A4.A4q.A1il"="#0099BB" ,"A4.A4ct.A4cs.A1il" = "#00AFC4","A1.A4.A1ac.A1bn"="#00C5CD")
ann_syms=colorRampPalette(c("#00AFC4","darkblue"))(length(unique(test$sym)))
#ann_syms=colorRampPalette(c("yellow","red"))(length(unique(test$sym)))
#ann_syms=colorRampPalette(c("pink","purple"))(length(unique(test$sym)))
names(ann_syms)<-unique(test$sym)[order(unique(test$sym))]
ann_colors=list(site=c(SI='#FF6663',BP='#FEB144',LK='#FDFD97',TK='#9EE09E',MI='#9EC1CF',OK='#CC99C9',CC7='black'))
ggplot(test)+
  geom_segment(aes(x=Axis.1,xend=sym.Axis.1,y=Axis.2,yend=sym.Axis.2,color=sym))+
  geom_point(aes(x=Axis.1,y=Axis.2,color=site),size=2)+
  scale_color_manual(values=c(ann_syms,ann_colors$site))+
  geom_point(aes(x=sym.Axis.1,sym.Axis.2,color=sym,size=2),shape=17)+
  theme_bw()
#facet_wrap(~sym)

# tanglegram
library(ape)
library(dendextend)
library(phytools)
# Convert trees to dendrograms
host_dend <- as.dendrogram(hclust(as.dist(z_final$H)))
para_dend <- as.dendrogram(hclust(as.dist(z_final$P)))
assoc <- z_final$HP@.Data
host.tree<-as.phylo(host_dend)
sym.tree<-as.phylo(para_dend)
cophyloplot(host.tree, sym.tree, assoc = assoc,show.tip.label = F)
library(Rtapas)
LFc <- max_cong(assoc[host.tree$tip.label,sym.tree$tip.label], host.tree, sym.tree, n=5, N=100, method = "paco", symmetric = TRUE,
                ei.correct = "sqrt.D", percentile = 0.01, res.fq = TRUE)
x11()
tangle_gram(host.tree,sym.tree,HS=assoc[host.tree$tip.label,sym.tree$tip.label],LFc,
            method="paco",colscale = "sequential", colgrad = c("darkblue","darkblue"),
            node.tag = F,link.lwd=1.5)
#dev.print(pdf,"/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/A_min5inds_Tanglegram_Jul2025.pdf") 
dev.off()

# plot out links
links<-as.data.frame(z_final$jackknife)
links$site=sapply(strsplit(rownames(links),"_"),FUN='[',1)
links$site=factor(links$site,levels=c("OK","MI","TK","LK","BP","SI"))
links$sym=gsub("^[0-9]+-","",gsub("OK-|MI-|TK-|LK-|BP-|SI-","",rownames(links)))

ggplot(links,aes(x=rownames(links),y=z_final$jackknife,fill=sym))+
  geom_bar(stat="identity")+
  scale_fill_manual(values=ann_syms)+
  facet_wrap(~sym)

ggplot(links,aes(x=site,y=z_final$jackknife))+
  geom_boxplot()+
  facet_wrap(~sym)

# also note MLL189host is messed up so need to redo that

