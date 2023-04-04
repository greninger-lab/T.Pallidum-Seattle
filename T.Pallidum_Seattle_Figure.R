library(tidyverse)
library(treeio)
library(ggtree)
library(phytools)
library(cowplot)
library(viridis)
library(phylotools)
library(Biostrings)
library(ggheatmap)
library(reshape2)
library(wesanderson)
library(ghibli)
library(dutchmasters)
library(gplots)
library(ggplot2)
#tree <- treeio::read.tree("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix/Masking_Proteins/Gubbins/Whole_Genome/masked_mafft_gubbins_WG.fasta.treefile")
#tree <- treeio::read.tree("/Users/administrator/Desktop/Matt_G_TP/Matt_G_SeattleSix_New_03032023/Masking_Regions/masked_regions/wholeGeome/WG_Tree/masked_mafft_gubbins_WG.fasta.treefile")

#tree <- treeio::read.tree("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/masking/masked_phylo/gubbins/tree/masked_mafft_gubbins_WG.fasta.treefile")

#no gubbins
#tree <- treeio::read.tree("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/masking/masked_phylo/gubbins/IQ_Tree_NoGubbins/combined.fasta.treefile")

#just deletion
#tree <- treeio::read.tree("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/masking/masked_phylo/gubbins/Just_Deletion/combined.fasta.treefile")

#No_Gub_WithSix
#tree <- treeio::read.tree("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/masking/masked_phylo/28Samp_Six_Seattle/combined_masked.fasta.treefile")

#Gub_JustSNP
#tree <- treeio::read.tree("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/masking/masked_phylo/gubbins/SNP_Gubbins_Tree/my_prefix.filtered_polymorphic_sites.fasta.treefile")

#NEwFixed Masked
#tree <- treeio::read.tree("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/Remask_Fix/Masked/GUB/masked_mafft_gubbins_WG.fasta.treefile")

#NEwFixed Masked_manual
tree <- treeio::read.tree("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/RemaskFix_Manual-Adjust_3_27_23/masked_mafft_gubbins_WG.fasta.treefile")

tree <- midpoint.root(tree)

tree$edge.length <- tree$edge.length * (1143784-17353)

#tree$edge.length <- tree$edge.length * (1143784-7407)

#iq tree no gubbins
#tree$edge.length <- tree$edge.length * (1143784-7052)

#just deletion
#tree$edge.length <- tree$edge.length * (1140399-27205)

fort <- fortify(tree)

#for labeling 95% bootstrap support
q <- ggtree(tree)
d <- q$data %>% filter(isTip == FALSE)# %>% arrange(desc(y))
#d <- d[!d$isTip,]

#d$label <- as.numeric(d$label)

d <- d[d$label > 95,]

x <- filter(fort, isTip == TRUE)
x <- arrange(x, y)
tips <- x[,c(4,7)]
q

ggtree_object <- ggtree(tree, size = 2.8) + theme_tree2() + geom_tiplab(linesize = 15,size =20) + labs(x="SNPs",size=55,face="bold") +
  #xlim(c(0,88)) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size=55,face="bold"),
    axis.title.x = element_text(size=55,face="bold"),
    axis.line =element_line(size=1.5)
  ) +
  geom_nodepoint(color="black", size = 4) +
  geom_tippoint(shape = "triangle", color = "grey", size = 4) 

#ggtree_object %<+% layout(rotate = list(node = 51, angle = 90))
#ggtree_object <- ggtree::rotate(ggtree_object, 56)

#print(ggtree_object)

#OGtips <- as.data.frame(tree$tip.label)

print(tips$label)

#metadata <- read.csv("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix/ggtree/metadata.csv",header=TRUE)
metadata <- read.csv("/Users/administrator/Desktop/Matt_G_TP/Matt_G_SeattleSix_New_03032023/Masking_Regions/masked_regions/wholeGeome/WG_Tree/ggtree/metadata.csv",header=TRUE)

#nogub+6
#metadata <- read.csv("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/masking/masked_phylo/28Samp_Six_Seattle/metadata.csv",header=TRUE)

# Create a vector of Sample names in the desired order
#sample_order <- c("P2220006", "P2120135", "P2220227", "R2110111", "P2120119", "R2110099", "UW15993L","P2220168", "P2120167", "UW15953", "UW15977O", "NC_021490","R2210139","P2220198", "V2140001", "UW15970L", "V2140014", "P2220282", "V2140005", "R2210209", "R2110052", "UW15939", "R2110091", "UW16057L", "P2120120", "P2120116", "NC_021508")
#sample_order <- c("UW15953", "P2120167", "P2220168", "R2110122", "UW15993L", "P2220227", "UW15977O", "P2220006", "P2120135", "R2110099", "R2110111", "P2120119", "P2120045", "R2210139", "P2220198", "NC_021490", "V2140014", "UW15970L", "V2140001", "V2140015", "P2220282", "V2140005", "R2210209", "R2110052", "UW15939", "R2110091","P2120120", "P2120116", "UW16057L", "NC_021508")

sample_order <- tips$label

#sample_order <- rev(sample_order)
# Convert the Sample variable to a factor with the desired order
metadata$Sample <- factor(metadata$Sample, levels = sample_order)

ggplot_object <- ggplot(metadata, aes(y = Sample, x = .5, fill = strain)) + 
  geom_tile(width = 0.5) +
  scale_fill_manual(values = c(
    `SS14` = "#B50A2A",
    `Nichols` = "#0E84B4")) +
  theme_minimal() +
  theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        #legend.position = "none",
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.title.x = element_text(angle = 90, vjust = .5, hjust=1,size=55,face="bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        legend.key.size = unit(1.5, unit="cm"),
        legend.title = element_text(size=60),
        legend.text = element_text(size=55)) +
  xlab("                            Lineage") +
 # theme(axis.title.x = element_text(hjust = 1))+
  guides(fill=guide_legend(title="Lineage"))

ggplot_object

ggplot_object2 <- ggplot(metadata, aes(y = Sample, x = .25, fill = Azithromycin_Resistance)) + 
  geom_tile() +
  scale_fill_manual(values = c(
    `Yes` = "#583B2B",
    `No` = "#AD8152")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("Azithromycin Resistance")

ggplot_object3 <- ggplot(metadata, aes(y = Sample, x = .25, fill = Sex)) + 
  geom_tile() +
  #scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest1")) +
  scale_fill_manual(values = c( `Unknown` = "white", `M` = "#E9D097", `F` = "#C5A387")) +

  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("                                  Sex") +
  guides(fill=guide_legend(title="Sex"))

ggplot_object4 <- ggplot(metadata, aes(y = Sample, x = .25, fill = MSM)) + 
  geom_tile() +

  theme_minimal() +
  #scale_fill_manual(values=wes_palette(n=3, name="GrandBudapest2")) +
  scale_fill_manual(values = c( `Unknown` = "white", `N` = "#1D2645", `Y` = "#B4DAE5", `N/A` = "#403369")) +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("                                 MSM") +
  guides(fill=guide_legend(title="MSM"))

ggplot_object4

ggplot_object5 <- ggplot(metadata, aes(y = Sample, x = .25, fill = tprk_donor_site_deletion)) + 
  geom_tile() +

  scale_fill_manual(values=wes_palette(n=2, name="Moonrise1")) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("Tprk Donor Site Deletion")+
  guides(fill=guide_legend(title="Tprk Donor Site Deletion"))

ggplot_object6 <- ggplot(metadata, aes(y = Sample, x = .25, fill = IV_drug_usage)) + 
  geom_tile() +
  #scale_fill_manual(values=wes_palette(n=3, name="Moonrise2")) +
  scale_fill_manual(values = c( `Unknown` = "white", `N` = "#26432F", `Y` = "#6FB382")) +

  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("                 IV Drug Usage")+
  guides(fill=guide_legend(title="IV Drug Usage"))

ggplot_object7 <- ggplot(metadata, aes(y = Sample, x = .25, fill = Meth_used_in_past_year)) + 
  geom_tile() +
  scale_fill_manual(values=wes_palette(n=3, name="Moonrise3")) +
  #scale_fill_manual(values = c( `Unknown` = "white")) +
 
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("Meth used in past year") +
  guides(fill=guide_legend(title="Meth used in past year"))

ggplot_object8 <- ggplot(metadata, aes(y = Sample, x = .25, fill = age_cat)) + 
  geom_tile() +
  scale_fill_manual(values = c(
    `Unknown` = "white",
    `20-24` = "gray20",
    `25-29` = "gray35",
    `30-34` = "gray50",
    `35-44` = "gray65",
    `45-54` = "gray80",
    `55-64` = "gray95")) +

  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("age cat")

ggplot_object9 <- ggplot(metadata, aes(y = Sample, x = .25, fill = race)) + 
  geom_tile() +
  #scale_fill_manual(values=wes_palette(n=5, name="BottleRocket1")) +
  scale_fill_dutchmasters(palette = "milkmaid") +

  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("Race")

ggplot_object10 <- ggplot(metadata, aes(y = Sample, x = .25, fill = Sample_Swab_Location)) + 
  geom_tile() +
  #scale_fill_manual(values=dutchmasters(n=9, name="little_street")) +
  #scale_fill_dutchmasters(palette = "pearl_earring") +
  
  scale_fill_manual(values = c( `Unknown` = "white", `Lesion` = "#06141F", `Throat` = "#742C14", `Rectal` = "#3D4F7D", `Vaginal` = "#CD4F38", `Papule near anus` = "#EAD890", `Oral` = "#3A160A")) +
  
  #scale_fill_manual(values = c( `Unknown` = "white")) +
  
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("   Sample Swab Location") +
  guides(fill=guide_legend(title="Sample Swab Location"))

print(ggplot_object10)

ggplot_object11 <- ggplot(metadata, aes(y = Sample, x = .25, fill = Disease_Stage)) + 
  geom_tile() +
  #scale_fill_dutchmasters(palette = "anatomy") +
  
  scale_fill_manual(values = c( `Unknown` = "white", `Late or unknown duration` = "grey", `Primary` = "#1F262E", `Early latent` = "#8F8093", `Secondary` = "#833437")) +
  
  
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    #legend.position = "none",
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(), 
    axis.title.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=55,face="bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.size = unit(1.5, unit="cm"),
    legend.title = element_text(size=60),
    legend.text = element_text(size=55)
  ) +
  xlab("                 Disease Stage") +
  guides(fill=guide_legend(title="Disease Stage"))

print(ggplot_object)
print(ggplot_object2)
print(ggplot_object3)
print(ggplot_object4)
print(ggplot_object5)
print(ggplot_object6)
print(ggplot_object7)
print(ggplot_object8)
print(ggplot_object9)
print(ggplot_object10)
print(ggplot_object11)
#fig2all <- plot_grid(ggtree_object, ggplot_object, ggplot_object2, ggplot_object3, ggplot_object4, ggplot_object5, ggplot_object6, ggplot_object7, ggplot_object8, ggplot_object9, ncol = 10, align = "h",rel_widths = c(1, .1, .1, .1, .1, .1, .1, .1, .1, .1))

p1 <-  get_legend(ggplot_object)
p2 <-  get_legend(ggplot_object2)
p3 <-  get_legend(ggplot_object3)
p4 <-  get_legend(ggplot_object4)
p5 <-  get_legend(ggplot_object5)
p6 <-  get_legend(ggplot_object6)
p7 <-  get_legend(ggplot_object7)
p8 <-  get_legend(ggplot_object8)
p9 <-  get_legend(ggplot_object9)
p10 <-  get_legend(ggplot_object10)
p11 <-  get_legend(ggplot_object11)

ggplot_object  <-  ggplot_object + theme(legend.position = "none")
ggplot_object2 <-  ggplot_object2 + theme(legend.position = "none")
ggplot_object3 <-  ggplot_object3 + theme(legend.position = "none")
ggplot_object4 <-  ggplot_object4 + theme(legend.position = "none")
ggplot_object5 <-  ggplot_object5 + theme(legend.position = "none")
ggplot_object6 <-  ggplot_object6 + theme(legend.position = "none")
ggplot_object7 <-  ggplot_object7 + theme(legend.position = "none")
ggplot_object8 <-  ggplot_object8 + theme(legend.position = "none")
ggplot_object9 <-  ggplot_object9 + theme(legend.position = "none")
ggplot_object10 <-  ggplot_object10 + theme(legend.position = "none")
ggplot_object11 <-  ggplot_object11 + theme(legend.position = "none")

#1 2 5
#4 6 7
#3 8 9

legendGroup <- plot_grid(p1,p2, p10, p5, ncol = 1, align = "v")
legendGroup2 <- plot_grid(p4, p3, p6, p11, ncol = 1, align = "v")

#fig2all <- plot_grid(ggtree_object, ggplot_object, ggplot_object2, ggplot_object10, ggplot_object5, legendGroup, ggplot_object4, ggplot_object3, ggplot_object6, ggplot_object11, legendGroup2, ncol = 11, align = "h",rel_widths = c(1, .02, .02, .02, .02, .07, .02, .02,.02, .02, .08))
fig2all <- plot_grid(ggtree_object, ggplot_object, ggplot_object2, ggplot_object10, ggplot_object5, ggplot_object4, ggplot_object3, ggplot_object6, ggplot_object11, ncol = 9, align = "h",rel_widths = c(1.3, .01, .01, .01, .01, .01, .01,.01, .01))

legendGroup3 <- plot_grid(p1, p2, p10, p5, p4, p3, p6, p11, ncol = 8, align = "h")
fig2all <- plot_grid(fig2all, legendGroup3, ncol = 1, align = "v", rel_heights = c(1, .07))

fig2all

#png(filename="~/Desktop/FIXED_MASKED_28+6.png", width=100, height=39, units="in", res=300)
png(filename="~/Desktop/FIXED_MASKED_28+6_3_27_23.png", width=100, height=100, units="in", res=300)
fig2all
dev.off()

SNP_Dist <- read.csv("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/RemaskFix_Manual-Adjust_3_27_23/SNPDist_molten.tsv", header = FALSE, sep = "\t", skip = 1)

SNP_Dist[, 1:2] <- t(apply(SNP_Dist[, 1:2], 1, function(x) sort(x)))

print(metadata)

# Remove duplicate rows
SNP_Dist <- unique(SNP_Dist)

print(SNP_Dist)

#merged_data <- merge(SNP_Dist, metadata, by = "Sample")

# Subset the data into two groups based on the strain column
ss14_data <- subset(metadata, strain == "SS14")
nicols_data <- subset(metadata, strain == "Nichols")

print(ss14_data)

SNP_Dist_NoReference <- SNP_Dist[!(SNP_Dist$V1 == "NC_021490" | SNP_Dist$V2 == "NC_021490"), ]
SNP_Dist_NoReference <- SNP_Dist_NoReference[!(SNP_Dist_NoReference$V1 == "NC_021508" | SNP_Dist_NoReference$V2 == "NC_021508"), ]

# filter the dataset to include only rows where both sample 1 and sample 2 are present in the dataset
filtered_dataset_ss14 <- SNP_Dist_NoReference[SNP_Dist_NoReference$V1 %in% ss14_data$Sample & SNP_Dist_NoReference$V2 %in% ss14_data$Sample, ]
print(filtered_dataset_ss14)

filtered_dataset_nicols_data <- SNP_Dist_NoReference[SNP_Dist_NoReference$V1 %in% nicols_data$Sample & SNP_Dist_NoReference$V2 %in% nicols_data$Sample, ]
print(filtered_dataset_nicols_data)

average_ss14 <- mean(filtered_dataset_ss14$V3)
average_nicols <- mean(filtered_dataset_nicols_data$V3)

print(average_ss14)

print(average_nicols)

filtered_dataset_nicols_data_noOutlier <- filtered_dataset_nicols_data[!(filtered_dataset_nicols_data$V1 == "R2210139" | filtered_dataset_nicols_data$V2 == "R2210139"), ]
filtered_dataset_nicols_data_noOutlier <- filtered_dataset_nicols_data_noOutlier[!(filtered_dataset_nicols_data_noOutlier$V1 == "P2220198" | filtered_dataset_nicols_data_noOutlier$V2 == "P2220198"), ]

average_nicols_noOutlier <- mean(filtered_dataset_nicols_data_noOutlier$V3)

print(average_nicols_noOutlier)

#remove sample pairs
filtered_dataset_nicols_data_noOutlier <- filtered_dataset_nicols_data_noOutlier[!(filtered_dataset_nicols_data_noOutlier$V1 == "R2110111" | filtered_dataset_nicols_data_noOutlier$V2 == "R2110111"), ]
filtered_dataset_ss14 <- filtered_dataset_ss14[!(filtered_dataset_ss14$V1 == "P2220282" | filtered_dataset_ss14$V2 == "P2220282"), ]

average_ss14 <- mean(filtered_dataset_ss14$V3)

print(average_ss14)

average_nicols_noOutlier <- mean(filtered_dataset_nicols_data_noOutlier$V3)

print(average_nicols_noOutlier)

library("heatmaply")

SNP_Dist_Heatmap <- read.csv("/Users/administrator/Desktop/Matt_G_TP/MattG_SeattleSix_New_030723/RemaskFix_Manual-Adjust_3_27_23/SNPDist.tsv", header = FALSE, sep = "\t")
SNP_Dist_Heatmap <- as.data.frame(SNP_Dist_Heatmap)
rownames(SNP_Dist_Heatmap) <- SNP_Dist_Heatmap[, 1]
SNP_Dist_Heatmap[, 1] <- NULL
colnames(SNP_Dist_Heatmap) <- SNP_Dist_Heatmap[1,]
SNP_Dist_Heatmap = SNP_Dist_Heatmap[-1,]
#SNP_Dist_Heatmap_2 <- as.matrix(SNP_Dist_Heatmap)
SNP_Dist_Heatmap_2 <- as.matrix(SNP_Dist_Heatmap)

#SNP_Dist_Heatmap_2 <- as.numeric(SNP_Dist_Heatmap_2)

mode(SNP_Dist_Heatmap_2) = "numeric"
data.frame(SNP_Dist_Heatmap_2)

#heatmap(SNP_Dist_Heatmap_2, col = colorRampPalette(c("white", "red"))(100))

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(-1,0,length=100),  # for red
               seq(0.01,0.8,length=100),           # for yellow
               seq(0.81,1,length=100))             # for green

my_palette <- colorRampPalette(c("white", "red", "darkslateblue"))(n = 299)
heatmap.2(SNP_Dist_Heatmap_2, notecol="black", cellnote = SNP_Dist_Heatmap_2,dendrogram = "none",key = FALSE, trace="none", col=my_palette)
#Colv="NA", breaks=col_breaks

#ggp <- ggplot(data) 
#ggp     
