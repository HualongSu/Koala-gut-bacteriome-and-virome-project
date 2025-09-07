

#
rm(list=ls())
##
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
###
setwd("")
# 
tree <- read.tree("input/version2.txt") # 
plot(tree)
# 
dt1 <- read.csv("input/taxonomy.csv")  # 
dt2 <- read.csv("input/ring1.csv")  # 
dt3 <- read.csv("input/ring2.csv")  # 
dt4 <- read.csv("input/barplot.csv")  # 
##
dt1$Phyla <- factor(dt1$Phyla, levels=c("Actinomycetota","Bacillota","Bacteroidota","Campylobacterota",
                                        "Cyanobacteriota","Deferribacterota","Desulfobacterota","Fusobacteriota",
                                        "Methanobacteriota","Planctomycetota","Pseudomonadota",
                                        "Synergistota","Verrucomicrobiota"))
# 
dt3$Type2 <- factor(dt3$Type2, levels=c("young", "adult",
                                        "middle-aged", "senior"))

##
# 
nodelab <- tree$node.label[nchar(tree$node.label)>0]  # 
nodeids <- nodeid(tree, nodelab)  # 

############################4、绘制环形树图##############

p1 <- ggtree(
  tree,
  layout="fan", open.angle = 90, # 
  size=0.15            # 
) 
p1
####################  5、添加门标签#############
# 
#colors <- c(
  #"#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", 
  #"#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", 
  #"#FC8D62", "#8DA0CB", "#E78AC3"
#)

colors <- c("#8C276A","#D7060F","#89B5E2","#C0D785","#EE7C70","#F9B85E",
            "#55A1B4","#1B5DA9","#9DD3C8","#7E6EB0","#E31880","#2E9A33")




p2 <- p1 %<+% dt1 +
  geom_tippoint(
    mapping=aes(
      fill=Phyla  #
    ),
    shape=21,  # 
    size=1.2,  #
    stroke=0.01,  # 
    position="identity",  # 
    show.legend=F  # 
  ) +
  scale_fill_manual(values=colors)  # 
p2
p2 <- p2 + new_scale_fill()  # 
p2
######################7、###############
p4 <- p2 +
  geom_fruit(
    data=dt3,
    geom=geom_tile,  # 
    mapping=aes(
      y=ID,  # 
      alpha=Abundance,  
      x=Type2,  # 
      fill=Type2  # 
    ),
    offset=0.0001,  # 
    pwidth=0.18   # 
  ) +
  scale_fill_manual(
    name="Fatty Acid metabolism",  # 
    values=c("#b22222", "#005500", "#0000be", "#9f1f9f"),  # 设置颜色
    guide=guide_legend(keywidth=0.35, keyheight=0.35, order=2)
  ) +
  scale_alpha_continuous(
    range=c(0, 1),  # 
    guide=guide_legend(keywidth=0.35, keyheight=0.35, order=3)
  ) +
  new_scale_fill()  # 
p4
##第二个版本
# 创建对数转换列
dt3$Abundance_log <- log10(dt3$Abundance + 1)  # 

# 作图
p4 <- p2 +
  geom_fruit(
    data = dt3,
    geom = geom_tile,  # 
    mapping = aes(
      y = ID,  
      alpha = Abundance_log,  
      x = Type2,  
      fill = Type2  
    ),
    offset = 0.0001,
    pwidth = 0.18
  ) +
  scale_fill_manual(
    name = "Group",
    values = rep("#064A93", length(unique(dt3$Type2))), 
    guide = guide_legend(keywidth = 0.35, keyheight = 0.35, order = 2)
  ) +
  scale_alpha_continuous(
    name = "log10(Abundance + 1)", 
    range = c(0.2, 1),  
    guide = guide_legend(keywidth = 0.35, keyheight = 0.35, order = 3)
  ) +
  new_scale_fill()

p4 <- p2 +
  geom_fruit(
    data = dt3,
    geom = geom_tile,  
    mapping = aes(
      y = ID,  
      alpha = Abundance,  
      x = Type2,  
      fill = Type2  
    ),
    offset = 0.0001,  
    pwidth = 0.18   
  ) +
  scale_fill_manual(
    name = "Fatty Acid metabolism",  
    values = rep("#064A93", length(unique(dt3$Type2))),  
    guide = guide_legend(keywidth = 0.35, keyheight = 0.35, order = 2)
  ) +
  scale_alpha_continuous(
    range = c(0, 1), 
    guide = guide_legend(keywidth = 0.35, keyheight = 0.35, order = 3)
  ) +
  new_scale_fill()  

p4
#######################8、最外层添加图形属性柱状图###########

p5 <- p4 +
  geom_fruit(data=dt4,
             geom=geom_bar,
             mapping=aes(
               y=ID,
               x=Length,
               fill=Phyla  
             ),
             stat="identity",
             orientation="y",  
             pwidth=0.5,
             position=position_dodgex()) +
  scale_fill_manual(
    values=colors,  
    guide=guide_legend(keywidth=0.35, keyheight=0.35, order=4)
  ) +
  geom_treescale(fontsize=1.2, linesize=0.3) +  
  theme(legend.position=c(0.93, 0.76),  
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"))  
p5
# 
ggsave("input/final2.pdf", dpi=800, plot=p5, width=10, height=10)  
ggsave("input/final2.tif", dpi=800, plot=p5, width=10, height=10) 
ggsave("input/final3.pdf", dpi=800, plot=p5, width=16, height=16)
ggsave("input/final4.pdf", dpi=800, plot=p5, width=10, height=1


