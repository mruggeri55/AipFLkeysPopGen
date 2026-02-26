setwd('/Users/maria/Desktop/Kenkel_lab/Aip_FLkeys/computational/2bRAD/2024/pop_gen/NGSadmix/final/')

all=readRDS('dfs/AipClean89_k2.Rdata')
noSIs=readRDS('dfs/bams50_noSIs_k3.Rdata')
SIsOnly=readRDS('dfs/SIsOnlyNoClones_k2.Rdata')

library(ggplot2)
library(forcats)

p.all=ggplot(all, aes(fill=admix, y=frac, x=ind)) + 
  geom_bar(position="stack", stat="identity",width=1) +
  ylab('ancestry')+
  xlab(NULL)+
  scale_fill_manual(values=c("#40C9A2","#F5B700"))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
  theme_bw(base_size = 18) + theme(legend.position = "none") + theme(axis.text.x=element_blank()) +
  facet_grid(~fct_inorder(pop), scales = "free", switch = "x", space = "free") +
  theme(panel.spacing.x = grid:::unit(0, "lines"), 
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.ticks.x = element_blank())

p.noSIs=ggplot(noSIs, aes(fill=admix, y=frac, x=ind)) + 
  geom_bar(position="stack", stat="identity",width=1) +
  ylab('ancestry')+
  xlab(NULL)+
  scale_fill_manual(values=c("#F5B700","#f9d466","#fce9b3"))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
  theme_bw(base_size = 18) + theme(legend.position = "none") + theme(axis.text.x=element_blank()) +
  facet_grid(~fct_inorder(pop), scales = "free", switch = "x", space = "free") +
  theme(panel.spacing.x = grid:::unit(0, "lines"), 
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.ticks.x = element_blank())

p.SIs=ggplot(SIsOnly, aes(fill=admix, y=frac, x=ind)) + 
  geom_bar(position="stack", stat="identity",width=1) +
  ylab('ancestry')+
  xlab(NULL)+
  scale_fill_manual(values=c("#40C9A2","#b5ffea"))+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(-.01,1.01), clip = "off")+
  theme_bw(base_size = 18) + theme(legend.position = "none") + theme(axis.text.x=element_blank()) +
  facet_grid(~fct_inorder(pop), scales = "free", switch = "x", space = "free") +
  theme(panel.spacing.x = grid:::unit(0, "lines"), 
        panel.border = element_rect(fill = NA, color = "black", linewidth = 1, linetype = "solid"),
        strip.background = element_rect(color = NA, fill = NA),
        axis.ticks.x = element_blank())

library(gridExtra)
bottom=grid.arrange(p.noSIs+labs(tag="B")+
                      theme(plot.tag=element_text(size=20),
                            plot.margin=unit(c(0,0.25,0.25,0.25),"cm")),
             p.SIs+labs(tag="C")+
               theme(axis.text.y = element_blank(),
                         axis.ticks.y=element_blank(),
                         axis.title.y = element_blank(),
                     plot.tag=element_text(size=20),
                     plot.margin=unit(c(0,0.25,0.25,0.25),"cm")),
             nrow=1,widths=c(1.5,1))

grid.arrange(p.all+labs(tag="A")+theme(plot.tag=element_text(size=20),
                                       plot.margin=unit(c(0.25,0.25,0,0.25),"cm")),
             bottom,nrow=2)




