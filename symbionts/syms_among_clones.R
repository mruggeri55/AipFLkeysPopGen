setwd('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025')

##### check for sym shuffling/switching
# read in sym relative abundance
zooxCounts = read.csv('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/its2_type_profiles/relAbundfromPhyloSeq.csv', row.names = 1)
zooxCounts$sample=rownames(zooxCounts)
# read in sample to MLL table
# using most conservative threshold for this analysis (aka individuals must be highly related to be considered clones)
#MLLs184=read.csv('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/MLL_lists/MLL184tab_pt023.csv',row.names = 1)
#colnames(MLLs184)[1]="name"
#MLLs184$name=gsub("_host.bam","",MLLs184$name)
# read in metadata
#meta=read.csv('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/clones/poppr_June2025/meta.csv',row.names = 1)
# merge meta with MLLs
#meta_MLL=merge(meta,MLLs184,by="name",all.y=T,all.x=T)
# export and fill in weird missing values
#write.csv(meta_MLL,"meta_MLLs184.csv")
meta_MLL=read.csv("syms_among_clones/meta_MLLs184.csv")
meta_MLL=meta_MLL[!grepl("TR2",meta_MLL$name),]
meta_MLL$MLL=as.factor(meta_MLL$MLL)

zooxCounts_master=merge(zooxCounts,meta_MLL,by='sample',all.x=T)
zooxCounts_master$sample[duplicated(zooxCounts_master$sample)]

rownames(zooxCounts_master)=zooxCounts_master$sample
zooxCounts_master2=zooxCounts_master[order(zooxCounts_master$MLL),]
zooxCounts_master2$order = c(1:nrow(zooxCounts_master2))

# note, only 272 samples with both high quality host and sym sequencing

library(tidyr)
library(dplyr)
melt_counts <- zooxCounts_master2 %>% gather(key=Reference,value=frac,colnames(zooxCounts_master2[2:36]))
melt_counts$Reference=as.factor(melt_counts$Reference)
levels(melt_counts$Reference)
melt_counts$site=as.character(melt_counts$site)
melt_counts$MLL=as.factor(melt_counts$MLL)
melt_counts$order=as.factor(melt_counts$order)

# count how many individuals are dominated (>50% abund) by each profile
sumDom = melt_counts[melt_counts$frac>=0.5,] %>% group_by(Reference) %>% summarise(n=n())

unique(melt_counts$sample[is.na(melt_counts$MLL)])
# some samples are missing MLLs do to failed preps/bad genotyping, exclude those

melt_counts$type = sapply(strsplit(as.character(melt_counts$Reference), split='\\.'),FUN='[',1)
melt_counts$genus = gsub("[[:digit:]]","",melt_counts$type)

# how many total clonal groups (aka where we could detect switching/shuffling)
temp=melt_counts[!is.na(melt_counts$MLL),] %>% group_by(MLL) %>% summarise(n=n_distinct(sample))
length(temp$MLL[temp$n>1])
clonal_groups=temp$MLL[temp$n>1]
# 42 clonal groups with more than one individual

# do dominant profs differ across ramets? types? genera?
DomSymsPerGenet = melt_counts[melt_counts$frac>=0.5 & !is.na(melt_counts$MLL),] %>% group_by(MLL) %>% summarise(Nprofs=n_distinct(Reference),Ntypes=n_distinct(type),Ngenera=n_distinct(genus))

# count how many MLLs are flexible at the profile, type, and genus level
length(DomSymsPerGenet$MLL[DomSymsPerGenet$Nprofs>1])
# 23 flexible genets (profile level) with new strict cutoff (184 MLLs)
length(DomSymsPerGenet$MLL[DomSymsPerGenet$Ntypes>1])
# 8 MLLs with different dominant profiles
length(DomSymsPerGenet$MLL[DomSymsPerGenet$Ngenera>1])
# 7 MLLs with different dominant genera

# make lists of flexible MLLs
SuperFlexGenets = DomSymsPerGenet$MLL[DomSymsPerGenet$Ngenera>1]
SemiFlexGenets = DomSymsPerGenet$MLL[DomSymsPerGenet$Ntypes>1]
# note only one MLL that differs at the type level (the rest are across genera)
FlexGenets = DomSymsPerGenet$MLL[DomSymsPerGenet$Nprofs>1]

#### plotting sym flexibility ####
library(ggplot2)
# A1 = colorRampPalette(c('dark green','navy'))
# A4 = colorRampPalette(c("navy", "turquoise3"))
# A = c(A1(4)[1:2],A4(11))
# B1 = colorRampPalette(c("red4", "#B2480D"))
# B19 = "gold"
# B2 = colorRampPalette(c("#D07E17", "#EEB422"))
# B = c(B1(6),B19,B2(4))
# C1='violet'
# C15 = colorRampPalette(c('darkorchid4','mediumorchid4'))
# C3 = colorRampPalette(c('hotpink','pink'))
# C = c(C1,C15(4),C3(4))
# D = c('dark grey','grey')
# colPalZoox=c(A,B,C,D)
colPalZoox=readRDS("/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/MS/final_figs/colPalZoox.RDS")

# Anew=colorRampPalette(c("#00AFC4","darkblue"))(length(unique(melt_counts$Reference[melt_counts$genus=="A"])))
# Bnew=colorRampPalette(c("red4","gold"))(length(unique(melt_counts$Reference[melt_counts$genus=="B"])))
# Cnew=colorRampPalette(c('pink','darkorchid4'))(length(unique(melt_counts$Reference[melt_counts$genus=="C"])))
# Dnew=colorRampPalette(c("dark grey","grey"))(length(unique(melt_counts$Reference[melt_counts$genus=="D"])))
# names(Anew)<-unique(melt_counts$Reference[melt_counts$genus=="A"])[order(unique(melt_counts$Reference[melt_counts$genus=="A"]))]
# names(Bnew)<-unique(melt_counts$Reference[melt_counts$genus=="B"])[order(unique(melt_counts$Reference[melt_counts$genus=="B"]))]
# names(Cnew)<-unique(melt_counts$Reference[melt_counts$genus=="C"])[order(unique(melt_counts$Reference[melt_counts$genus=="C"]))]
# names(Dnew)<-unique(melt_counts$Reference[melt_counts$genus=="D"])[order(unique(melt_counts$Reference[melt_counts$genus=="D"]))]
# colPalZoox=c(Anew,Bnew,Cnew,Dnew)

plot_list = list()
for (i in SuperFlexGenets) {
  p=ggplot(na.omit(melt_counts[melt_counts$MLL==i,]), aes(fill=Reference, y=frac, x=name)) + 
    geom_bar(position="stack", stat="identity") +
    labs(subtitle=paste('MLL',i) ,y="")+
    xlab(NULL)+
    scale_fill_manual(values=colPalZoox)+
    coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
    theme_bw(base_size = 16) + 
    theme(axis.text.x=element_blank(),legend.position = 'none',axis.ticks.x = element_blank()) +
    facet_grid(~sub_rt, scales = "free", switch = "x", space = "free")+
    theme(strip.text = element_text(size = 8,angle=90))
  
  plot_list[[i]] = p
}
library(cowplot)
plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
          plot_list[[5]],plot_list[[6]],plot_list[[7]],
          ncol=3)

# plot all together but separated by bolding
library(forcats)
ggplot(na.omit(melt_counts[melt_counts$MLL %in% SuperFlexGenets,]), aes(fill=Reference, y=frac, x=name)) + 
  geom_bar(position="stack", stat="identity",width=1) +
  ylab('relative abundance')+
  xlab('clonal group')+
  scale_fill_manual(values=colPalZoox)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
  theme_bw(base_size = 18) + 
  theme(legend.position = "none") + 
  theme(axis.text.x=element_blank()) +
  facet_grid(~fct_inorder(MLL), scales = "free", switch = "x", space = "free") +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.ticks.x = element_blank())

# plot all semi flex genets
plot_list = list()
for (i in c(SemiFlexGenets[!SemiFlexGenets %in% SuperFlexGenets],FlexGenets[!FlexGenets %in% SuperFlexGenets])) {
  p=ggplot(na.omit(melt_counts[melt_counts$MLL==i,]), aes(fill=Reference, y=frac, x=name)) + 
    geom_bar(position="stack", stat="identity") +
    labs(subtitle=paste('MLL',i) ,y="")+
    xlab(NULL)+
    scale_fill_manual(values=colPalZoox)+
    coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
    theme_bw(base_size = 16) + theme(axis.text.x=element_blank(),legend.position = 'none')
  plot_list[[i]] = p
}

plot_grid(plot_list[[1]],plot_list[[2]],plot_list[[3]],plot_list[[4]],
          plot_list[[5]],plot_list[[6]],plot_list[[7]],plot_list[[8]],plot_list[[9]],
          ncol=3)

plot_grid(plot_list[[10]],plot_list[[11]],plot_list[[12]],
          plot_list[[13]],plot_list[[14]],plot_list[[15]],
          plot_list[[16]],
          ncol=3)

#plot_grid(plot_list[[19]],plot_list[[20]],plot_list[[21]],
#          plot_list[[22]],plot_list[[23]],
#          ncol=3)

ggplot(na.omit(melt_counts[melt_counts$MLL %in% c(SemiFlexGenets[!SemiFlexGenets %in% SuperFlexGenets],FlexGenets[!FlexGenets %in% SuperFlexGenets]),]), 
       aes(fill=Reference, y=frac, x=name)) + 
  geom_bar(position="stack", stat="identity",width=1) +
  ylab('relative abundance')+
  xlab('clonal group')+
  scale_fill_manual(values=colPalZoox)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
  theme_bw(base_size = 18) + 
  theme(legend.position = "none") + 
  theme(axis.text.x=element_blank()) +
  facet_grid(~fct_inorder(MLL), scales = "free", switch = "x", space = "free") +
  theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.ticks.x = element_blank())

# plot for all clonal groups
plot_list = list()
site=c("OK","MI","TK","LK","BP","SI")
scale=c()
for (i in site) {
  p=ggplot(na.omit(melt_counts[melt_counts$MLL %in% clonal_groups & melt_counts$site==i,]), 
       aes(fill=Reference, y=frac, x=name)) + 
  geom_bar(position="stack", stat="identity",width=1) +
  labs(y = i,x = NULL) +
  scale_fill_manual(values=colPalZoox)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),breaks=c(0,0.5,1)) +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
  theme_bw(base_size = 18) + 
  theme(legend.position = "none") + 
  facet_grid(~fct_inorder(MLL), scales = "free", switch = "x", space = "free") +
  theme(panel.spacing.x = grid:::unit(0.5, "lines"), 
        panel.border = element_rect(fill = NA, color = "black", size = 0.5, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.text.x=element_blank(),
        strip.text.x = element_blank())
  plot_list[[i]] = p
  scale=c(scale,
          (length(unique(melt_counts$name[melt_counts$MLL %in% clonal_groups & melt_counts$site==i]))+
             length(unique(melt_counts$MLL[melt_counts$MLL %in% clonal_groups & melt_counts$site==i]))-1))
}
plot_grid(plot_grid(plot_list[[1]],NULL,rel_widths=c(scale[[1]],46-scale[[1]])),
          plot_grid(plot_list[[2]],NULL,rel_widths=c(scale[[2]],40-scale[[2]])),
          plot_grid(plot_list[[3]],NULL,rel_widths=c(scale[[3]],46-scale[[3]])),
          plot_grid(plot_list[[4]],NULL,rel_widths=c(scale[[4]],46-scale[[4]])),
          plot_grid(plot_list[[5]],NULL,rel_widths=c(scale[[5]],46-scale[[5]])),
          plot_grid(plot_list[[6]],NULL,rel_widths=c(scale[[6]],35-scale[[6]])),
          ncol=1)
# this is close but not perfect
