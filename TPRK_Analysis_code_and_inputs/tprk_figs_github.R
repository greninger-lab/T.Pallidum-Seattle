library(tidyverse)
library(phylotools)
library(reshape2)
library(cowplot)
library(viridis)
library(paletteer)
library(treeio)
library(ggtree)
library(phytools)

`%!in%` <- negate(`%in%`)

seqs <- read_csv("~/Desktop/FINAL_avg_data_seqs.csv")
#colnames(seqs) <- gsub("_extA", "", colnames(seqs))
#colnames(seqs) <- gsub("_extB", "", colnames(seqs))

bin <- as.data.frame(ifelse(seqs[,4:ncol(seqs)] > 0, 1, 0))
seqs_per_sample <- data.frame(Sample=colnames(bin), unique_seqs=colSums(bin))
rownames(seqs_per_sample) <- NULL

meta <- read_csv("~/Desktop/per_sample_meta.csv") 
bin_meta <- inner_join(meta, seqs_per_sample)


cors <- read_csv("/Users/Nicole/Desktop/per_sample_tprk_pearsons.csv")

#Supp figures

cor.test(bin_meta$input_copies, bin_meta$unique_seqs) #p-value = 0.2969, cor 0.2043703 

suppA <- ggplot(bin_meta, aes(x=input_copies, y=unique_seqs)) +
  geom_point(size=3, shape=21, fill="grey", color="black") + 
  xlim(c(0,2500)) +
  ylim(c(0,80))+
  xlab(expression(paste("Input ",italic("tp47 "), "copies"))) +
  ylab("Unique V Sequences") +
  theme_classic()
suppA

x <- bin_meta %>% filter(Lineage == "Nichols") #`Disease Stage` == "Secondary")            #
y <- bin_meta %>% filter(Lineage == "SS14")  #`Disease Stage` == "Primary")            
t.test(x$unique_seqs, y$unique_seqs)
#p for Nichols vs SS14: 0.5859
#p for primary secondary: 0.2827

#Non parametric: 
wilcox.test(x$unique_seqs, y$unique_seqs)
#p for primary secondary: 0.5235
#p for Nichols vs SS14: 0.9815
suppB <- ggplot(bin_meta %>% filter(`Disease Stage` == "Primary" | `Disease Stage` == "Secondary"), aes(x=`Disease Stage`, y=unique_seqs, fill=Lineage)) +
  geom_boxplot() + 
  ylim(c(0,80))+
  ylab("Unique V Sequences") +
  scale_fill_manual(values=c("#0E84B4", "#B50A2A")) +
  theme_classic()
suppB

suppAB <- plot_grid(suppA, suppB, labels=c("A", "B"), rel_widths=c(1,1.3))
suppAB

pdf("~/Desktop/tprk_supp.pdf", width=6.5, height=3.25)
suppAB
dev.off()



#Fig 3
#first need to add in 
tree <- treeio::read.tree("/Users/Nicole/Desktop/seattle_tree.newick")
tree <- midpoint.root(tree)
tree$edge.length <- tree$edge.length * (1143784-28602)
fort <- fortify(tree)

q <- ggtree(tree)
d <- q$data %>% filter(isTip == FALSE)# %>% arrange(desc(y))
x <- filter(fort, isTip == TRUE)
x <- arrange(x, y)
tips <- x[,c(4,7)]
q
d$label <- as.numeric(d$label)
d <- d[d$label > 90,]

tree <- drop.tip(tree, c("NC_021508", "NC_021490"))

ggtree_object <- ggtree(tree, size = .5) + 
  theme_tree2() +
  #geom_tiplab(align = TRUE, linesize=.25, size = 3) +
 # geom_tiplab(align = TRUE, linesize=.25, size = 3) +
  geom_point2(aes(label=label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 90),size=1, color = "red") +
#  geom_text2(aes(subset = label == "R-22-10139", color = "red")) +
  labs(x="SNPs",size=6,face="bold") +
  coord_cartesian(clip="off")+
  xlim(c(0,260)) +
  #scale_x_continuous(breaks = seq(0, 250, by = 50)) +
  theme_void() +
  theme(
   
    plot.margin = margin(0, 1.5, 0, 0, "cm")
  ) 
print(ggtree_object)

print(d$node)

ggtree_object <- ggtree::rotate(ggtree_object, 40)

print(tips$label)

pdf(file = "~/Desktop/tree.pdf", width = 1.5, height = 3.25)
ggtree_object
dev.off()


#then modified slightly in illustrator


cors_melt <- melt(cors)
colnames(cors_melt) <- c("SampleA", "SampleB", "Cor. Coeff")
cors_melt$SampleA <- factor(cors_melt$SampleA, levels=meta$Sample)
cors_melt$SampleB <- factor(cors_melt$SampleB, levels=rev(meta$Sample))

A <- ggplot(cors_melt, aes(x=SampleA, y=SampleB, fill=`Cor. Coeff`)) +
  geom_tile() +
  scale_fill_viridis(option="C") +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
A
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 5, spaceLegend = 0.3) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

# Apply on original plot
addSmallLegend(A)
A


#now add part B: 

cor_stats <- read_csv("/Users/Nicole/Desktop/corr_long_form.csv")
cor_stats$type <- str_to_title(cor_stats$type)
cor_stats <- cor_stats %>% filter(cor_stats$SampleA != cor_stats$SampleB)

#all the stats
type_aov <- aov(coeff ~ type, data=cor_stats, contrasts=TRUE)
summary(type_aov) #2e-16
tukey.test <- TukeyHSD(type_aov)
tukey.test

cor_stats$type <- factor(cor_stats$type, levels=c("Same Patient", "Linked Cluster", "Same Lineage", "Opposite Lineage"))

plot3B <- ggplot(cor_stats, aes(x=type, y=coeff, fill=type)) +
  geom_boxplot(outlier.shape=NA, color="black") +
  geom_jitter(width=0.2, size=2, alpha=0.5, shape=21, color="black") +
  ylab(expression(paste(italic("tprK "), "Correlation Coefficient"))) +
  geom_segment(x=1, xend=1.95, y=1.02, yend=1.02) +
  geom_text(aes(label = "****", x=1.5, y=1.05)) +
  geom_segment(x=2.05, xend=2.95, y=1.02, yend=1.02) +
  geom_text(aes(label = "****", x=2.5, y=1.05)) +
  geom_segment(x=3.05, xend=4, y=1.02, yend=1.02) +
  geom_text(aes(label = "*", x=3.5, y=1.05)) +
  scale_fill_viridis_d(option="C", direction=-1) +
  scale_x_discrete(labels=str_wrap(levels(cor_stats$type),width=9)) +
  theme_classic(base_size=7) +
  theme(legend.position="none",
        axis.title.x=element_blank())
plot3B


fig3AB <- plot_grid(NULL, A, plot3B, ncol=3, rel_widths=c(0.75,3.5,2.25), labels=c("A", "", "B"))

pdf(file="~/Desktop/fig3AB.pdf", width=5.5, height=3.25)
fig3AB
dev.off()




#fig 5: 

dsko <- bin_meta %>% filter(Lineage == "Nichols" & `Disease Stage` == "Secondary")
dsko_x <- dsko %>% filter(`tprK Donor Site Deletion` == "Yes")
dsko_y <- dsko %>% filter(`tprK Donor Site Deletion` == "No")
t.test(dsko_x$unique_seqs, dsko_y$unique_seqs) #p-value = 0.006041 mean of x mean of y 17.00000  34.57143 
#may be diffs based on randomization in downsampling

fig5A <- ggplot(bin_meta %>% filter(Lineage == "Nichols" & `Disease Stage` == "Secondary"), aes(x=`tprK Donor Site Deletion`, y=unique_seqs)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, size=3, alpha=0.5, ) +
  xlab(expression(paste(italic("tprK "), "Donor Site Deletion"))) +
  ylab("Unique V Sequences") +
  geom_segment(x=1, xend=2, y=max(dsko$unique_seqs)+5, yend=max(dsko$unique_seqs)+5) +
  geom_text(aes(label = "*", x=1.5, y=max(dsko$unique_seqs)+6)) +
  ylim(c(0,65)) +
  theme_classic(base_size=7) 
fig5A

dsko_use_n <- read_csv("~/Desktop/dsko_use_n.csv")

ds_no <- dsko_use_n %>% filter(`tprK Donor Site Deletion`== "No")
ds_yes <- dsko_use_n %>% filter(`tprK Donor Site Deletion`== "Yes")
xxx <- t.test(ds_yes$KO_percent, ds_no$KO_percent) #p-value = 0.01521 mean of x mean of y  5.47619  13.90821 
xxx

fig5b <- ggplot(dsko_use_n, aes(x=`tprK Donor Site Deletion`, y=KO_percent)) + #`tprK Donor Site Deletion`
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.2, height =0, size=3, alpha=0.5) +
  xlab(expression(paste(italic("tprK "), "Donor Site Deletion"))) +
  ylab("% of Unique V Sequences Using \nDS Within Deletion") +
  geom_segment(x=1, xend=2, y=max(dsko_use_n$KO_percent)+5, yend=max(dsko_use_n$KO_percent)+5) +
  geom_text(aes(label = "*", x=1.5, y=max(dsko_use_n$KO_percent)+6)) +
  ylim(c(0,40)) +
  theme_classic(base_size=7)
fig5b

ds <- read_csv("/Users/Nicole/Desktop/donor_site_usage.csv")
ds_meta <- ds %>% select(Name, Region, DS, Minimum)
ds_melt <- melt(ds, id.vars = c("Name", "Region", "DS", "Minimum"))
ds_melt$Name <- factor(ds_melt$Name, levels=unique(ds_melt$Name))
ds_melt$variable <- factor(ds_melt$variable, levels=rev(meta$Sample))

ds_melt <- ds_melt %>% filter(value > 0)

fig5c <- ggplot(ds_melt, aes(x=Name, y=variable, fill=value)) +
  geom_tile(color="white", fill="darkseagreen1") +
  geom_text(aes(label=value), size=3) +
  geom_segment(aes(x=28.5, xend=36.5, y=24.5, yend=24.5), linetype="dashed", color="red") +
  geom_segment(aes(x=28.5, xend=36.5, y=28.9, yend=28.9), linetype="dashed", color="red") +
  
  geom_segment(aes(x=28.5, xend=28.5, y=24.5, yend=29), linetype="dashed", color="red") +
  geom_segment(aes(x=36.5, xend=36.5, y=24.5, yend=29), linetype="dashed", color="red") +
  xlab("Donor Site (arranged by genomic position)") +
  scale_x_discrete(position = "top") +
  theme_classic(base_size=8) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle=45, hjust=0, vjust=),
        axis.title.y = element_blank())
fig5c

Cspace <- plot_grid(fig5c, NULL, ncol=2, rel_widths = c(40,1), labels="C")
Cspace

ds_sum <- data.frame(Name=ds$Name, Total=rowSums(ds[,2:29]))
ds_sum$Name <- factor(ds_sum$Name, levels=ds_sum$Name)
which_not_used <- ds_sum %>% filter(Total == 0)
ds_sum <- ds_sum %>% filter(Total >0)
pt2_5c <- ggplot(ds_sum, aes(x=Name, y=1, fill=Total)) +
  geom_tile(color="white") +
  geom_text(aes(label=Total), size=2.5) +
  scale_fill_paletteer_c("grDevices::Light Grays", direction = -1) +
  theme_void() +
  theme(legend.position="none")
pt2_5c

C2space <- plot_grid(NULL, pt2_5c, NULL, ncol=3, rel_widths = c(0.95,10,0.4))
C2space


fig5C_total <- plot_grid(Cspace, C2space, ncol=1, rel_heights=c(10,0.75))
fig5C_total

fig5AB <- plot_grid(fig5A, fig5b, ncol=2, labels=c("A","B"))

fig5 <- plot_grid(fig5AB, fig5C_total, ncol=1)
fig5


pdf(file="~/Desktop/Fig5.pdf", width=6.5, height=6.5)
fig5
dev.off()
