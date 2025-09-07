rm(list=ls())#clear Global Environment
setwd("")
##
library(tidyverse)
da <- read.table("output_count_koala_bins.txt", header = T, sep = "\t", row.names = 1)
data_group = read.delim("group.txt", row.names = 1)
da %>%  head()

#
otu <- da %>%

  as.data.frame()

# otu
#
#install.packages("microeco")
library(microeco)
data <- microtable$new(otu_table = otu,
                       sample_table = data_group)
data$otu_table
##
# 1.4 
#install.packages("remotes")    # install.packages("devtools")
#remotes::install_github("ChiLiubio/mecodev")
library(mecodev)
colSums(data$otu_table) 
max(colSums(data$otu_table)) 
min(colSums(data$otu_table)) 

set.seed(125)

otu_nom <- as.data.frame(
  t(rrarefy(t(otu), min(colSums(otu))))
) %>%
  filter(., rowSums(.) > 0)


# 
min(colSums(otu_nom))
#
##

write.csv(
  otu_nom,
  file = "otu_filt.csv",
  row.names = FALSE
)
data_group = read.delim("group.txt", row.names = 1)

data <- microtable$new(otu_table = otu_nom,                       
                       sample_table = data_group)
data$otu_table
#data$cal_alphadiv(measures = NULL, PD = FALSE)
a1 <- trans_alpha$new(dataset = data)
alpha <- a1$data_alpha %>%  
  pivot_wider(names_from = Measure,              
              values_from = Value)

group_counts <- table(alpha$group)

valid_groups <- names(group_counts[group_counts >= 2])


alpha <- alpha[alpha$group %in% valid_groups, ]
library(ggplot2)
library(ggsignif)
col <- c("#da0000","#6600cc","#3D9F3C","#9ED17B","#ff9900","#6a3d00","#9DC7DD","#367DB0")
scale_color_manual(values = col)


groups <- c("Young", "Adult", "Senior")


comparisons <- combn(groups, 2, simplify = FALSE)

library(ggplot2)
library(dplyr)
#install.packages("rstatix")
library(rstatix)
#install.packages("multcompView")
library(multcompView)
library(ggplot2)
library(dplyr)
library(agricolae)


alpha_clean <- alpha %>% filter(!is.na(group), !is.na(Shannon))


anova_result_shannon <- aov(Shannon ~ group, data = alpha_clean)
summary(anova_result_shannon)


lsd_result_shannon <- LSD.test(anova_result_shannon, "group", p.adj = "bonferroni")


group_letters_shannon <- lsd_result_shannon$groups
group_letters_shannon$group <- rownames(group_letters_shannon)
colnames(group_letters_shannon)[which(names(group_letters_shannon) == "groups")] <- "label"


label_df_shannon <- alpha_clean %>%
  group_by(group) %>%
  summarise(y = max(Shannon, na.rm = TRUE) + 0.2) %>%
  left_join(group_letters_shannon[, c("group", "label")], by = "group")
alpha_clean$group <- factor(alpha_clean$group, levels = c("Young", "Adult", "Senior"))

max(alpha_clean$Shannon)
min(alpha_clean$Shannon)

p1 <- ggplot(alpha_clean, aes(x = group, y = Shannon, color = group)) +
  geom_violin(fill = NA, scale = "width", linewidth = 0.8, alpha = 0.7) +
  geom_jitter(size = 1, alpha = 0.7, width = 0.2) +
  stat_summary(fun = median, geom = "point", shape = 3, size = 3, colour = "black") +
  geom_text(data = label_df_shannon, aes(x = group, y = y, label = label),
            color = "black", size = 5) +
  scale_color_manual(values = col) +
  scale_y_continuous(limits = c(5, 8), breaks = seq(5, 8, 0.5)) +
  labs(x = NULL, y = "Shannon diversity", color = NULL) +
  theme_bw() +
  theme(
    axis.text         = element_text(size = 9, color = "black"),
    axis.title.y      = element_text(size = 9, angle = 90),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", size = 1),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.title        = element_text(hjust = 0.5),
    legend.position   = "none"
  )

p1





p1.1 <- p1 + annotate("text", x = 1, y = max(alpha_clean$Shannon) + 0.3, 
              label = "P = 0.001", hjust = 0, size = 5)
p1.1


anova_result_InvSimpson <- aov(InvSimpson ~ group, data = alpha_clean)
summary(anova_result_InvSimpson)

lsd_result_InvSimpson <- LSD.test(anova_result_InvSimpson, "group", p.adj = "bonferroni")


group_letters_InvSimpson <- lsd_result_InvSimpson$groups
group_letters_InvSimpson$group <- rownames(group_letters_InvSimpson)
colnames(group_letters_InvSimpson)[which(names(group_letters_InvSimpson) == "groups")] <- "label"


label_df_InvSimpson <- alpha_clean %>%
  group_by(group) %>%
  summarise(y = max(InvSimpson, na.rm = TRUE) + 25) %>%
  left_join(group_letters_InvSimpson[, c("group", "label")], by = "group")


max(alpha_clean$InvSimpson)
min(alpha_clean$InvSimpson)
p2 <- ggplot(alpha_clean, aes(x = group, y = InvSimpson, color = group)) +
  geom_violin(fill = NA, scale = "width", linewidth = 0.8, alpha = 0.7) +
  geom_jitter(size = 1, alpha = 0.7, width = 0.2) +
  stat_summary(fun = median, geom = "point", shape = 3, size = 3, colour = "black") +
  geom_text(data = label_df_InvSimpson, aes(x = group, y = y, label = label),
            color = "black", size = 5) +
  scale_color_manual(values = col) +
  scale_y_continuous(limits = c(80, 700), breaks = seq(80, 700, 50)) +
  labs(x = NULL, y = "InvSimpson diversity", color = NULL) +
  theme_bw() +
  theme(
    axis.text         = element_text(size = 9, color = "black"),
    axis.title.y      = element_text(size = 9, angle = 90),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", size = 1),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.title        = element_text(hjust = 0.5),
    legend.position   = "none"
  )

p2
p2.1 <- p2 + annotate("text", x = 1, y = max(alpha_clean$InvSimpson) + 8, 
              label = "P = 0.001", hjust = 0, size = 5)
p2.1



anova_result_Simpson <- aov(Simpson ~ group, data = alpha_clean)
summary(anova_result_Simpson)


lsd_result_Simpson <- LSD.test(anova_result_Simpson, "group", p.adj = "bonferroni")


group_letters_Simpson <- lsd_result_Simpson$groups
group_letters_Simpson$group <- rownames(group_letters_Simpson)
colnames(group_letters_Simpson)[which(names(group_letters_Simpson) == "groups")] <- "label"


label_df_Simpson <- alpha_clean %>%
  group_by(group) %>%
  summarise(y = max(Simpson, na.rm = TRUE) + 0.05) %>%
  left_join(group_letters_Simpson[, c("group", "label")], by = "group")



max(alpha_clean$Simpson)
min(alpha_clean$Simpson)
p3 <- ggplot(alpha_clean, aes(x = group, y = Simpson, color = group)) +
  geom_violin(fill = NA, scale = "width", linewidth = 0.8, alpha = 0.7) +
  geom_jitter(size = 1, alpha = 0.7, width = 0.2) +
  stat_summary(fun = median, geom = "point", shape = 3, size = 3, colour = "black") +
  geom_text(data = label_df_Simpson, aes(x = group, y = y, label = label),
            color = "black", size = 5) +
  scale_color_manual(values = col) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.1)) +
  labs(x = NULL, y = "Simpson diversity", color = NULL) +
  theme_bw() +
  theme(
    axis.text         = element_text(size = 9, color = "black"),
    axis.title.y      = element_text(size = 9, angle = 90),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", size = 1),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.title        = element_text(hjust = 0.5),
    legend.position   = "none"
  )

p3


max(alpha_clean$Pielou)
min(alpha_clean$Pielou)

anova_result_Pielou <- aov(Pielou ~ group, data = alpha_clean)
summary(anova_result_Pielou)


lsd_result_Pielou <- LSD.test(anova_result_Pielou, "group", p.adj = "bonferroni")


group_letters_Pielou <- lsd_result_Pielou$groups
group_letters_Pielou$group <- rownames(group_letters_Pielou)
colnames(group_letters_Pielou)[which(names(group_letters_Pielou) == "groups")] <- "label"


label_df_Pielou <- alpha_clean %>%
  group_by(group) %>%
  summarise(y = max(Pielou, na.rm = TRUE) + 0.03) %>%  
  left_join(group_letters_Pielou[, c("group", "label")], by = "group")


p4 <- ggplot(alpha_clean, aes(x = group, y = Pielou, color = group)) +
  geom_violin(fill = NA, scale = "width", linewidth = 0.8, alpha = 0.7) +
  geom_jitter(size = 1, alpha = 0.7, width = 0.2) +
  stat_summary(fun = median, geom = "point", shape = 3, size = 3, colour = "black") +
  geom_text(data = label_df_Pielou, aes(x = group, y = y, label = label),
            color = "black", size = 5) +
  scale_color_manual(values = col) +
  scale_y_continuous(
    limits = c(0.5, 0.8), 
    breaks = seq(0.5, 0.8, 0.1)
  ) +
  labs(x = NULL, y = "Pielou diversity", color = NULL) +
  theme_bw() +
  theme(
    axis.text         = element_text(size = 9, color = "black"),
    axis.title.y      = element_text(size = 9, angle = 90),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", size = 1),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.title        = element_text(hjust = 0.5),
    legend.position   = "none"
  )


p4

p4.1 <- p4 + annotate("text", x = 1, y = max(alpha_clean$Pielou) + 0.1, 
                      label = "P = 0.001", hjust = 0, size = 5)
p4.1
##Fisher

max(alpha_clean$Fisher)
min(alpha_clean$Fisher)
# Step 1: 
anova_result_Fisher <- aov(Fisher ~ group, data = alpha_clean)
summary(anova_result_Fisher)

# Step 2:
lsd_result_Fisher <- LSD.test(anova_result_Fisher, "group", p.adj = "bonferroni")

# Step 3: 
group_letters_Fisher <- lsd_result_Fisher$groups
group_letters_Fisher$group <- rownames(group_letters_Fisher)
colnames(group_letters_Fisher)[which(names(group_letters_Fisher) == "groups")] <- "label"

# Step 4: 
label_df_Fisher <- alpha_clean %>%
  group_by(group) %>%
  summarise(y = max(Fisher, na.rm = TRUE) + 200) %>%  
  left_join(group_letters_Fisher[, c("group", "label")], by = "group")

# Step 5: 
p5 <- ggplot(alpha_clean, aes(x = group, y = Fisher, color = group)) +
  geom_violin(fill = NA, scale = "width", linewidth = 0.8, alpha = 0.7) +
  geom_jitter(size = 1, alpha = 0.7, width = 0.2) +
  stat_summary(fun = median, geom = "point", shape = 3, size = 3, colour = "black") +
  geom_text(data = label_df_Fisher, aes(x = group, y = y, label = label),
            color = "black", size = 5) +
  scale_color_manual(values = col) +
  scale_y_continuous(
    limits = c(2000,6500), 
    breaks = seq(2000, 65000, 500)
  ) +
  labs(x = NULL, y = "Fisher diversity", color = NULL) +
  theme_bw() +
  theme(
    axis.text         = element_text(size = 9, color = "black"),
    axis.title.y      = element_text(size = 9, angle = 90),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", size = 1),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.title        = element_text(hjust = 0.5),
    legend.position   = "none"
  )


p5
p5.1 <- p5 + annotate("text", x = 1, y = max(alpha_clean$Fisher) + 0.1, 
              label = "P = 0.001", hjust = 0, size = 5)
p5.1

##Chao1
max(alpha_clean$Chao1)
min(alpha_clean$Chao1)


# Step 1:
anova_result_Chao1 <- aov(Chao1 ~ group, data = alpha_clean)
summary(anova_result_Chao1)

# Step 2: 
lsd_result_Chao1 <- LSD.test(anova_result_Chao1, "group", p.adj = "bonferroni")

# Step 3:
group_letters_Chao1 <- lsd_result_Chao1$groups
group_letters_Chao1$group <- rownames(group_letters_Chao1)
colnames(group_letters_Chao1)[which(names(group_letters_Chao1) == "groups")] <- "label"

# Step 4: 
label_df_Chao1 <- alpha_clean %>%
  group_by(group) %>%
  summarise(y = max(Chao1, na.rm = TRUE) + 700) %>%  
  left_join(group_letters_Chao1[, c("group", "label")], by = "group")

# Step 5: 
p6 <- ggplot(alpha_clean, aes(x = group, y = Chao1, color = group)) +
  geom_violin(fill = NA, scale = "width", linewidth = 0.8, alpha = 0.7) +
  geom_jitter(size = 1, alpha = 0.7, width = 0.2) +
  stat_summary(fun = median, geom = "point", shape = 3, size = 3, colour = "black") +
  geom_text(data = label_df_Chao1, aes(x = group, y = y, label = label),
            color = "black", size = 5) +
  scale_color_manual(values = col) +
  scale_y_continuous(
    limits = c(21000,55000), 
    breaks = seq(21000, 55000, 4000)
  ) +
  labs(x = NULL, y = "Chao1 diversity", color = NULL) +
  theme_bw() +
  theme(
    axis.text         = element_text(size = 9, color = "black"),
    axis.title.y      = element_text(size = 9, angle = 90),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", size = 1),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.title        = element_text(hjust = 0.5),
    legend.position   = "none"
  )

# 
p6

p6.1 <- p6 + annotate("text", x = 1, y = max(alpha_clean$Chao1) + 0.1, 
              label = "P = 0.001", hjust = 0, size = 5)

p6.1
##
###Richness画图，也就是observed那一列的值
max(alpha_clean$Observed)
min(alpha_clean$Observed)


# Step 1: 
anova_result_Observed <- aov(Observed ~ group, data = alpha_clean)
summary(anova_result_Observed)

# Step 2: 
lsd_result_Observed <- LSD.test(anova_result_Observed, "group", p.adj = "bonferroni")

# Step 3:
group_letters_Observed <- lsd_result_Observed$groups
group_letters_Observed$group <- rownames(group_letters_Observed)
colnames(group_letters_Observed)[which(names(group_letters_Observed) == "groups")] <- "label"

# Step 4: 
label_df_Observed <- alpha_clean %>%
  group_by(group) %>%
  summarise(y = max(Observed, na.rm = TRUE) + 1000) %>%  
  left_join(group_letters_Observed[, c("group", "label")], by = "group")

# Step 5: 
p7 <- ggplot(alpha_clean, aes(x = group, y = Observed, color = group)) +
  geom_violin(fill = NA, scale = "width", linewidth = 0.8, alpha = 0.7) +
  geom_jitter(size = 1, alpha = 0.7, width = 0.2) +
  stat_summary(fun = median, geom = "point", shape = 3, size = 3, colour = "black") +
  geom_text(data = label_df_Observed, aes(x = group, y = y, label = label),
            color = "black", size = 5) +
  scale_color_manual(values = col) +
  scale_y_continuous(
    limits = c(19000,54000),  # 修正此处
    breaks = seq(19000,54000, 5000)
  ) +
  labs(x = NULL, y = "Richness", color = NULL) +
  theme_bw() +
  theme(
    axis.text         = element_text(size = 9, color = "black"),
    axis.title.y      = element_text(size = 9, angle = 90),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", size = 1),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.title        = element_text(hjust = 0.5),
    legend.position   = "none"
  )


p7
p7.1 <- p7 + annotate("text", x = 1, y = max(alpha_clean$Observed) + 15, 
                      label = "P = 0.001", hjust = 0, size = 5)

p7.1
###计算Pielou
###计算Pielou
max(alpha_clean$Pielou)
min(alpha_clean$Pielou)


# Step 1: 
anova_result_Pielou <- aov(Pielou ~ group, data = alpha_clean)
summary(anova_result_Pielou)

# Step 2: 
lsd_result_Pielou <- LSD.test(anova_result_Pielou, "group", p.adj = "bonferroni")

# Step 3: 
group_letters_Pielou <- lsd_result_Pielou$groups
group_letters_Pielou$group <- rownames(group_letters_Pielou)
colnames(group_letters_Pielou)[which(names(group_letters_Pielou) == "groups")] <- "label"

# Step 4: 
label_df_Pielou <- alpha_clean %>%
  group_by(group) %>%
  summarise(y = max(Pielou, na.rm = TRUE) + 0.05) %>% 
  left_join(group_letters_Pielou[, c("group", "label")], by = "group")

# Step 5: 
p8 <- ggplot(alpha_clean, aes(x = group, y = Pielou, color = group)) +
  geom_violin(fill = NA, scale = "width", linewidth = 0.8, alpha = 0.7) +
  geom_jitter(size = 1, alpha = 0.7, width = 0.2) +
  stat_summary(fun = median, geom = "point", shape = 3, size = 3, colour = "black") +
  geom_text(data = label_df_Pielou, aes(x = group, y = y, label = label),
            color = "black", size = 5) +
  scale_color_manual(values = col) +
  scale_y_continuous(
    limits = c(0.5,0.8),  # 修正此处
    breaks = seq(0.5, 0.8, 0.05)
  ) +
  labs(x = NULL, y = "Pielou", color = NULL) +
  theme_bw() +
  theme(
    axis.text         = element_text(size = 9, color = "black"),
    axis.title.y      = element_text(size = 9, angle = 90),
    panel.grid        = element_blank(),
    panel.border      = element_rect(color = "black", size = 1),
    panel.background  = element_rect(fill = "white", colour = NA),
    plot.title        = element_text(hjust = 0.5),
    legend.position   = "none"
  )

# 输出图形
p8
p8.1 <- p8 + annotate("text", x = 1, y = max(alpha_clean$Pielou) + 0.05, 
                      label = "P = 0.001", hjust = 0, size = 5)

p8.1






###beta
#### library packages and load the mag and group
library(vroom)
library(ggplot2)
library(vegan)
library(tidyverse)
library(dplyr)
t_otu_nom <- t(otu_nom)
#### caculate distance matrix of all samples
set.seed(123)
t_otu_nom_bray <- vegdist(t_otu_nom,method = "bray")
t_otu_nom_bray

### anosim of samples of differnet trenches
anosim_result_dis <- anosim(t_otu_nom_bray, data_group$group, permutations = 999)
summary(anosim_result_dis) #get R and p of anosim

permanova_res2 <- adonis2(
  t_otu_nom_bray ~ group,
  data        = data_group,
  permutations = 999
)

print(permanova_res2)
###1. NMDS
###1.1 NMDS of Trench
nmds <- metaMDS(t_otu_nom_bray, k = 2)
stress <- nmds$stress ### stress
stress
nm_df <- as.data.frame(nmds$points)#get points
nm_df$sample <- row.names(nm_df)  
group<-read.table("group.txt",
                  sep='\t', header=T,check.names=FALSE )
nm_df1 <- merge(nm_df,group,by="sample") # all sample points

nm_df1 <- merge(nm_df,group,by="sample") # all sample points
p9<-ggplot(data=nm_df1,aes(MDS1, MDS2))+#
  theme_bw()+
  theme(panel.grid = element_blank()) +
  geom_point(aes(fill=group),size=2.5,color="NA",shape=21,alpha=1)+#
  labs(x=paste0("NMDS1 (Stress=0.150)"),
       y=paste0("NMDS2 "))+
  scale_fill_manual(values = col)+
  theme(axis.title.x=element_text(size=9,color="black"),#
        axis.title.y=element_text(size=9,angle=90,color="black"),#
        axis.text.y=element_text(size=9,color="black"),#x
        axis.text.x=element_text(size=9,color="black"),
        panel.border = element_rect(size = 1, color = "black"),
        legend.position   = "none"     
  )

print(p9)
p9.1 <- p9 + annotate("text", x = 0.5, y = 1.8, 
                      label = "P = 0.001", hjust = 0, size = 5)

p9.1
#####

pcoa <- cmdscale (t_otu_nom_bray,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
pc12$sample <- row.names(pc12)
head(pc12)
eig = summary(eigenvals(pcoa))#over all
axis = paste0("PCoA", 1:ncol(eig))
eig = data.frame(Axis = axis, t(eig)[, -3])
#get pco1 and pco2 of all sample

pco1 = round(eig[1, 3] * 100, 2)
pco2 = round(eig[2, 3] * 100, 2)
pco1
pco2
pcoa_all <- merge(pc12,group,by="sample")

###
adonis_result <- adonis2(t_otu_nom_bray ~ group$group, permutations = 999)
print(adonis_result)
##

p10<-ggplot(data=pcoa_all,aes(x=V1,y=V2))+#
  theme_bw()+
  theme(panel.grid = element_blank()) +
  geom_point(aes(fill=group),size=2.5,color="NA",shape=21,alpha=1)+
  labs(x=paste0("PCo1 = 29.65%"),
       y=paste0("PCo2 = 23.17%"))+
  scale_fill_manual(values = col)+
  theme(axis.title.x=element_text(size=9,color="black"),#
        axis.title.y=element_text(size=9,angle=90,color="black"),#
        axis.text.y=element_text(size=9,color="black"),#x
        axis.text.x=element_text(size=9,color="black"),
        panel.border = element_rect(size = 1, color = "black"),
        legend.position   = "none"     # 关闭右侧图例
  )
print(p10)

p10.1 <- p10 + annotate("text", x = -0.2, y = 0.3, 
                      label = "P = 0.001", hjust = 0, size = 5)
p10.1

#ggsave("p10.pdf", plot = p10, width = 10, height = 10, dpi = 1200)

##
library(patchwork)
groups_combined_alpha_beta_pd <- (p1.1 + p2.1 + p5.1) /
  (p6.1 + p7.1+ p8.1) / (p9.1 + p10.1)

groups_combined_alpha_beta_pd

ggsave("groups_alpha_beta_mirus.pdf", plot = groups_combined_alpha_beta_pd, width = 10, height = 10, dpi = 1200)

##
write.csv(
  alpha_clean,
  file = "alpha_clean.csv",
  row.names = FALSE
)

##
# 如果 group_age 是以 Sample 为 row.names 的 dataframe，先把行名变为一列
group_age = read.delim("group_age.txt", row.names = 1)
group_age$Sample <- rownames(group_age)
# 合并时用完整的 group_age 而不是 group_age$group
sandian_df <- merge(alpha_clean, group_age, by = "Sample")

######绘图#######

sandian_df$group <- factor(sandian_df$group, levels = c("Young", "Adult", "Senior"))
##age2Richness
##
library(ggpmisc)

p11 <- ggplot(sandian_df, aes(x = age, y = Observed)) +
  geom_point(aes(color = group), shape = 20, size = 3, alpha = 1) + # 修改 shape，移除边框
  geom_smooth(
    method = "lm", se = FALSE, color = "grey", 
    formula = y ~ x,
    linetype = "dashed", alpha = 0.5, linewidth = 0.8 # 修改为虚线
  ) +
  geom_smooth(
    method = "lm", se = FALSE, color = "black", 
    formula = y ~ poly(x, 2, raw = TRUE),
    alpha = 0.5, linewidth = 0.8 # 修改为点线
  ) +
  stat_poly_eq(
    formula = y ~ x, 
    aes(label = paste( 
      after_stat(rr.label), 
      after_stat(p.value.label),
      sep = "~~~")), 
    parse = TRUE, label.x = 0.95, label.y = 0.95, color = "grey"
  ) +
  stat_poly_eq(data = sandian_df,
               formula = y ~ poly(x, 2, raw = TRUE), 
               aes(label = paste(
                 after_stat(rr.label),
                 after_stat(p.value.label),
                 sep = "~~~")), 
               parse = TRUE, label.x = 0.95, label.y = 0.9,
               color = "black") +
  scale_color_manual(values = c("#da0000","#6600cc","#3D9F3C","#9ED17B","#ff9900","#6a3d00","#9DC7DD","#367DB0")) + # 修改点的颜色
  scale_x_continuous(breaks = seq(3, 10, by = 1),limits = c(3, 10))+##x轴的刻度问题
  scale_y_continuous(breaks = seq(19000,54000, by = 5000),limits = c(19000,54000))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA, color = "black", size = 0.8, linetype = "solid"),
        axis.text.x = element_text(color = "black"), # 设置 x 轴刻度标签的文本大小和颜色
        axis.text.y = element_text(color = "black"),# 设置 y 轴刻度标签的文本大小和颜色
        legend.position = "right") +  # 去掉图例
  labs(x = "Age", y = "Richness")
p11


##age2shannon


p12 <- ggplot(sandian_df, aes(x = age, y = Shannon)) +
  geom_point(aes(color = group), shape = 20, size = 3, alpha = 1) + # 修改 shape，移除边框
  geom_smooth(
    method = "lm", se = FALSE, color = "grey", 
    formula = y ~ x,
    linetype = "dashed", alpha = 0.5, linewidth = 0.8 # 修改为虚线
  ) +
  geom_smooth(
    method = "lm", se = FALSE, color = "black", 
    formula = y ~ poly(x, 2, raw = TRUE),
    alpha = 0.5, linewidth = 0.8 # 修改为点线
  ) +
  stat_poly_eq(
    formula = y ~ x, 
    aes(label = paste( 
      after_stat(rr.label), 
      after_stat(p.value.label),
      sep = "~~~")), 
    parse = TRUE, label.x = 0.95, label.y = 0.95, color = "grey"
  ) +
  stat_poly_eq(data = sandian_df,
               formula = y ~ poly(x, 2, raw = TRUE), 
               aes(label = paste(
                 after_stat(rr.label),
                 after_stat(p.value.label),
                 sep = "~~~")), 
               parse = TRUE, label.x = 0.95, label.y = 0.9,
               color = "black") +
  scale_color_manual(values = c("#da0000","#6600cc","#3D9F3C","#9ED17B","#ff9900","#6a3d00","#9DC7DD","#367DB0")) + # 修改点的颜色
  scale_x_continuous(breaks = seq(3, 10, by = 1),limits = c(3, 10))+##x轴的刻度问题
  scale_y_continuous(breaks = seq(5, 8, by = 0.5),limits = c(5, 8))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA, color = "black", size = 0.8, linetype = "solid"),
        axis.text.x = element_text(color = "black"), # 设置 x 轴刻度标签的文本大小和颜色
        axis.text.y = element_text(color = "black"),# 设置 y 轴刻度标签的文本大小和颜色
        legend.position = "right") +  # 去掉图例
  labs(x = "Age", y = "Shannon")
p12


ggsave("age2Shannon.pdf", plot = p12, width = 10, height = 10, dpi = 1200)


##age2InvSimpson

max(sandian_df$InvSimpson)
min(sandian_df$InvSimpson)
p13 <- ggplot(sandian_df, aes(x = age, y = InvSimpson)) +
  geom_point(aes(color = group), shape = 20, size = 3, alpha = 1) + # 修改 shape，移除边框
  geom_smooth(
    method = "lm", se = FALSE, color = "grey", 
    formula = y ~ x,
    linetype = "dashed", alpha = 0.5, linewidth = 0.8 # 修改为虚线
  ) +
  geom_smooth(
    method = "lm", se = FALSE, color = "black", 
    formula = y ~ poly(x, 2, raw = TRUE),
    alpha = 0.5, linewidth = 0.8 # 修改为点线
  ) +
  stat_poly_eq(
    formula = y ~ x, 
    aes(label = paste( 
      after_stat(rr.label), 
      after_stat(p.value.label),
      sep = "~~~")), 
    parse = TRUE, label.x = 0.95, label.y = 0.95, color = "grey"
  ) +
  stat_poly_eq(data = sandian_df,
               formula = y ~ poly(x, 2, raw = TRUE), 
               aes(label = paste(
                 after_stat(rr.label),
                 after_stat(p.value.label),
                 sep = "~~~")), 
               parse = TRUE, label.x = 0.95, label.y = 0.9,
               color = "black") +
  scale_color_manual(values = c("#da0000","#6600cc","#3D9F3C","#9ED17B","#ff9900","#6a3d00","#9DC7DD","#367DB0")) + # 修改点的颜色
  scale_x_continuous(breaks = seq(3, 10, by = 1),limits = c(3, 10))+##x轴的刻度问题
  scale_y_continuous(breaks = seq(80, 700, by = 50),limits = c(80, 700))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA, color = "black", size = 0.8, linetype = "solid"),
        axis.text.x = element_text(color = "black"), # 设置 x 轴刻度标签的文本大小和颜色
        axis.text.y = element_text(color = "black"),# 设置 y 轴刻度标签的文本大小和颜色
        legend.position = "right") +  # 去掉图例
  labs(x = "Age", y = "Invsimpson")
p13


###age2Pielou
max(alpha_clean$Pielou)
min(alpha_clean$Pielou)
p14 <- ggplot(sandian_df, aes(x = age, y = Pielou)) +
  geom_point(aes(color = group), shape = 20, size = 3, alpha = 1) + # 修改 shape，移除边框
  geom_smooth(
    method = "lm", se = FALSE, color = "grey", 
    formula = y ~ x,
    linetype = "dashed", alpha = 0.5, linewidth = 0.8 # 修改为虚线
  ) +
  geom_smooth(
    method = "lm", se = FALSE, color = "black", 
    formula = y ~ poly(x, 2, raw = TRUE),
    alpha = 0.5, linewidth = 0.8 # 修改为点线
  ) +
  stat_poly_eq(
    formula = y ~ x, 
    aes(label = paste( 
      after_stat(rr.label), 
      after_stat(p.value.label),
      sep = "~~~")), 
    parse = TRUE, label.x = 0.95, label.y = 0.95, color = "grey"
  ) +
  stat_poly_eq(data = sandian_df,
               formula = y ~ poly(x, 2, raw = TRUE), 
               aes(label = paste(
                 after_stat(rr.label),
                 after_stat(p.value.label),
                 sep = "~~~")), 
               parse = TRUE, label.x = 0.95, label.y = 0.9,
               color = "black") +
  scale_color_manual(values = c("#da0000","#6600cc","#3D9F3C","#9ED17B","#ff9900","#6a3d00","#9DC7DD","#367DB0")) + # 修改点的颜色
  scale_x_continuous(breaks = seq(3, 10, by = 1),limits = c(3, 10))+##x轴的刻度问题
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1),limits = c(0.5, 1))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(fill = NA, color = "black", size = 0.8, linetype = "solid"),
        axis.text.x = element_text(color = "black"), # 设置 x 轴刻度标签的文本大小和颜色
        axis.text.y = element_text(color = "black"),# 设置 y 轴刻度标签的文本大小和颜色
        legend.position = "right") +  # 去掉图例
  labs(x = "Age", y = "Pielou")
p14

##
merge_age2_diver <- (p12 + p13 + p14) 
 

merge_age2_diver

ggsave("merge_age2_diver.pdf", plot = merge_age2_diver, width = 10, height = 10, dpi = 1200)


