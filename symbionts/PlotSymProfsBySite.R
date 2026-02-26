# plot sym profs by site
setwd('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/AipFLkeysSymVar/CatRunsOut/its2_type_profiles')
# merge zooxcounts with metadata
zooxCounts = read.csv('relAbundfromPhyloSeq.csv', row.names = 1)
zooxCounts$sample=rownames(zooxCounts)
meta=read.csv('../meta_finalwGenets_nodups.csv',row.names = 1)

zooxCounts_master=merge(zooxCounts,meta,by='sample')
rownames(zooxCounts_master)=zooxCounts_master$sample
#zooxCounts_sub=zooxCounts_master[,c(1:36,38:42)]
#zooxCounts=zooxCounts_sub

library(RColorBrewer)
library(tidyr)

zooxCounts_master2=zooxCounts_master[order(zooxCounts_master$genotype),]
zooxCounts_master2$order = c(1:nrow(zooxCounts_master2))

melt_counts <- zooxCounts_master2 %>% gather(key=Reference,value=frac,colnames(zooxCounts_master2[2:36]))
melt_counts <- melt_counts[!melt_counts$site=="CC7",]

melt_counts$site=factor(melt_counts$site,levels=c('OK','MI','TK','LK','BP','SI'))
# reorder by subsite and root
melt_counts_order=melt_counts[order(melt_counts$site,melt_counts$subsite,melt_counts$root),]
melt_counts=melt_counts_order

melt_counts$Reference=as.factor(melt_counts$Reference)
levels(melt_counts$Reference)

A1 = colorRampPalette(c('dark green','navy'))
A4 = colorRampPalette(c("navy", "turquoise3"))
A = c(A1(4)[1:2],A4(11))
#Bfun = colorRampPalette(c("red4", 'goldenrod2'))
B1 = colorRampPalette(c("red4", "#B2480D"))
B19 = "gold"
B2 = colorRampPalette(c("#D07E17", "#EEB422"))
B = c(B1(6),B19,B2(4))
C1='violet'
C15 = colorRampPalette(c('darkorchid4','mediumorchid4'))
C3 = colorRampPalette(c('hotpink','pink'))
C = c(C1,C15(4),C3(4))
D = c('dark grey','grey')
colPalZoox=c(A,B,C,D)

melt_counts$genus=sapply(strsplit(as.character(melt_counts$Reference),"[0-9]"),FUN='[',1)
Anew=colorRampPalette(c("#00AFC4","darkblue"))(length(unique(melt_counts$Reference[melt_counts$genus=="A"])))
Bnew=colorRampPalette(c("red4","gold"))(length(unique(melt_counts$Reference[melt_counts$genus=="B"])))
Cnew=colorRampPalette(c('pink','darkorchid4'))(length(unique(melt_counts$Reference[melt_counts$genus=="C"])))
Dnew=colorRampPalette(c("dark grey","grey"))(length(unique(melt_counts$Reference[melt_counts$genus=="D"])))
names(Anew)<-unique(melt_counts$Reference[melt_counts$genus=="A"])[order(unique(melt_counts$Reference[melt_counts$genus=="A"]))]
names(Bnew)<-unique(melt_counts$Reference[melt_counts$genus=="B"])[order(unique(melt_counts$Reference[melt_counts$genus=="B"]))]
names(Cnew)<-unique(melt_counts$Reference[melt_counts$genus=="C"])[order(unique(melt_counts$Reference[melt_counts$genus=="C"]))]
names(Dnew)<-unique(melt_counts$Reference[melt_counts$genus=="D"])[order(unique(melt_counts$Reference[melt_counts$genus=="D"]))]
colPalZoox=c(Anew,Bnew,Cnew,Dnew)

# plot by site and subsite
plot_list = list()
for (i in levels(melt_counts$site)) {
  p=ggplot(melt_counts[melt_counts$site==i,], aes(fill=Reference, y=frac, x=name)) + 
    geom_bar(position="stack", stat="identity",width=1) +
    scale_fill_manual(values=colPalZoox)+
    labs(y = i,x = NULL) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0),breaks=c(0,0.5,1)) +
    coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
    theme_bw(base_size = 18) + 
    theme(legend.position = "none") + 
    theme(axis.text.x=element_blank()) +
    facet_grid(~fct_inorder(subsite), scales = "free", switch = "x", space = "free") +
    theme(panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
          strip.background = element_rect(color = NA, fill = NA),
          axis.ticks.x = element_blank())+
    theme(axis.title.y = element_text(angle=0,vjust=0.5))
  
  plot_list[[i]] = p
}

legend=get_legend(ggplot(melt_counts[melt_counts$site=="OK",], aes(fill=Reference, y=frac, x=name)) + 
                    geom_bar(position="stack", stat="identity",width=1) +
                    scale_fill_manual(values=colPalZoox))

library(ggpubr)
library(grid)
figure=ggarrange(ggarrange(plot_list[[1]],NULL,plot_list[[2]],NULL,plot_list[[3]],NULL,plot_list[[4]],NULL,plot_list[[5]],NULL,plot_list[[6]],NULL,
                           ncol=1,heights=c(rep(c(1,-0.15),times=11))),legend, ncol=2, widths = c(2, 1))
X11()
annotate_figure(figure, left = textGrob("relative abundance", rot = 90, vjust = 1, gp = gpar(cex = 1.5)))
dev.print(pdf,"/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/syms/SymProfsBySite_newcolors_13Aug2025.pdf") 
dev.off()

