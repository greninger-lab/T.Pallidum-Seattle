library(tidyverse)
library(Biostrings)
library(viridisLite)
library(vroom)
library(reshape2)

#setwd("~/Desktop/")
options(scipen = 5)
`%!in%` <- negate(`%in%`)

meta <- read_csv("~/Desktop/tprk_meta.csv") 

########################################
####TO GENERATE READS FROM RAW INPUT FILES:

path <- "/Users/Nicole/Desktop/count_data/" #need final / for vroom to work
setwd(path) #UPDATED
files <- list.files(path=path, pattern="_final_data.csv")
alldata <- vroom(files, id="fastq_name")
alldata$fastq_name <- gsub("_final_data.csv", "", alldata$fastq_name)

all_abs <- alldata[,c(1,2,3,5)]
all_abs <- all_abs %>% filter(fastq_name %in% meta$fastq_name)
length(unique(all_abs$fastq_name)) #should be same length as #in metadata
fastqs <- unique(all_abs$fastq_name)

reads <- all_abs
reads <- reads[!is.na(reads$Read),]

#remove reads with ambiguous bases
reads <- reads[!grepl("N", reads$Read),] 

names <- meta %>% select(fastq_name, SampleName)
reads <- inner_join(reads, names)
reads <- reads %>% select(-fastq_name)

cast <- dcast(reads, Region + Read ~ SampleName, value.var="Count") 
keep <- rowSums(cast[3:ncol(cast)], na.rm=TRUE) > 0
sum(keep) #expect same
cast <- cast[keep,]

reads <- cast

#translate
reads_info <- reads %>% select(Region, Read)
reads_info$AAseq <- as.character(translate(DNAStringSet(reads_info$Read)))
reads_info$total_counts <- rowSums(reads[3:ncol(reads)], na.rm=TRUE)
reads_info$V_length <- str_length(reads_info$Read)
#reads_info$InFrame <- ifelse(reads_info$V_length %% 3, FALSE, TRUE) #is V region divisible by 3? 
#This doesn't work because V4 region not correct length. For V4, need to subtract the last 2 AA (almost always GG, not considered part of V region....)
reads_info$V_length <- ifelse(reads_info$Region == "V4", reads_info$V_length - 2, reads_info$V_length) 
#now try again.... 
reads_info$InFrame <- ifelse(reads_info$V_length %% 3, FALSE, TRUE) #is V region divisible by 3? 


reads_info$Nonsense <- str_detect(reads_info$AAseq, "\\*")

valid_reads <- cbind(reads_info, reads[3:ncol(reads)])
valid_reads <- valid_reads %>% filter(InFrame == TRUE & Nonsense == FALSE & total_counts > 0)
valid_reads <- valid_reads[,-c(4:7)]

zeroed <- valid_reads[,4:ncol(valid_reads)] 
zeroed[is.na(zeroed)] <- 0
sum(rowSums(zeroed) > 0) #same length as zeroed
zeroed <- cbind(valid_reads[,1:3], zeroed) 

per_v <- zeroed %>% select(-c(Read, AAseq))
nv <- aggregate(. ~ Region, data=per_v, sum)

#how many reads to use?
nlines <- min(nv[2:ncol(nv)]) #3471

#assuming error rate of ~0.5%, is filtering #?
error_threshold <- round(nlines*0.005) #17. 


valid_reads2 <- valid_reads
valid_reads2$key <- paste0(valid_reads2$Region, "_", valid_reads2$Read, "_", valid_reads2$AAseq)
ds <- data.frame(key="x")

samples <- meta[meta$SampleName %in% colnames(nv), ]
#i <- "P_22_20198_B0403"
#set.seed(24) 
for (i in samples$SampleName) {
  print(i)
  a <- valid_reads2 %>% select(Region, key, one_of(i))
  a <- a[complete.cases(a),]
  av1 <- a %>% filter(Region == "V1") 
  longv1 <- av1[rep(row.names(av1), av1[[i]]), 1:ncol(av1)] 
  v1len <- min(nrow(longv1), nlines)
  dsav1 <- longv1[sample(nrow(longv1), v1len, replace=TRUE), ]
  av2 <- a %>% filter(Region == "V2")
  longv2 <- av2[rep(row.names(av2), av2[[i]]), 1:ncol(av2)]
  v2len <- min(nrow(longv2), nlines)
  dsav2 <- longv2[sample(nrow(longv2), v2len, replace=TRUE), ]
  av3 <- a %>% filter(Region == "V3")
  longv3 <- av3[rep(row.names(av3), av3[[i]]), 1:ncol(av3)]
  v3len <- min(nrow(longv3), nlines)
  dsav3 <- longv3[sample(nrow(longv3), v3len, replace=TRUE), ]
  av4 <- a %>% filter(Region == "V4")
  longv4 <- av4[rep(row.names(av4), av4[[i]]), 1:ncol(av4)]
  v4len <- min(nrow(longv4), nlines)
  dsav4 <- longv4[sample(nrow(longv4), v4len, replace=TRUE), ]
  av5 <- a %>% filter(Region == "V5")
  longv5 <- av5[rep(row.names(av5), av5[[i]]), 1:ncol(av5)]
  v5len <- min(nrow(longv5), nlines)
  dsav5 <- longv5[sample(nrow(longv5), v5len, replace=TRUE), ]
  av6 <- a %>% filter(Region == "V6")
  longv6 <- av6[rep(row.names(av6), av6[[i]]), 1:ncol(av6)]
  v6len <- min(nrow(longv6), nlines)
  dsav6 <- longv6[sample(nrow(longv6), v6len, replace=TRUE), ]
  av7 <- a %>% filter(Region == "V7")
  longv7 <- av7[rep(row.names(av7), av7[[i]]), 1:ncol(av7)]
  v7len <- min(nrow(longv7), nlines)
  dsav7 <- longv7[sample(nrow(longv7), v7len, replace=TRUE), ]
  subV <- rbind(dsav1, dsav2, dsav3, dsav4, dsav5, dsav6, dsav7)
  agg <- subV %>% group_by(key) %>% summarize(n=n()) 
  split <- as.data.frame(str_split_fixed(agg$key, "_", 3))
  agg <- cbind(split, agg)
  agg <- agg %>% arrange(V1, desc(n)) 
  agg <- agg %>% select(key, n)
  colnames(agg) <- c("key", i)
  ds <- full_join(ds, agg)
}


ds <- ds %>% filter(key !="x")
ds[is.na(ds)] <- 0
non_zero <- rowSums(ds[2:ncol(ds)]) > 0
ds <- ds[non_zero,]
total_reads_used <- as.data.frame(colSums(ds[2:ncol(ds)])) #all same, 7*nlines

ds_split <- as.data.frame(str_split_fixed(ds$key, "_", 3))
colnames(ds_split) <- c("Region", "Read", "AAseq")
ds <- cbind(ds_split, ds) %>% arrange(Region)

#write_csv(ds, "~/Desktop/tprk_replicates_downsample3471_final.csv")



####NOW: THESE ARE ALL VALID READS, use for figures 
repA <- samples %>% dplyr::filter(Replicate == "A")
repB <- samples %>% dplyr::filter(Replicate == "B")

identical(repA$SampleExt, repB$SampleExt)

repA <- repA %>% filter(repA$SampleExt %in% repB$SampleExt)
repB <- repB %>% filter(repB$SampleExt %in% repA$SampleExt)

identical(repA$SampleExt, repB$SampleExt) #should be true now

cors <- data.frame(SampleExt=repA$SampleExt, Pearson="")
#reads_sample_only <- reads[,4:ncol(reads)]

#i <- 17
for (i in 1:nrow(cors)) {
  x <- cors[i,1]
  xA <- paste0(x, "_A")
  xB <- paste0(x, "_B")
  y <- cor.test(ds[[xA]], ds[[xB]])
  z <- y$estimate
  cors[i,2] <- z
}


identical(cors$SampleExt, repA$SampleExt) #true
cors$Pearson <- as.numeric(cors$Pearson)


for_avg <- cors %>% filter(Pearson > 0.85) 

sample_order <- meta %>% filter(Replicate=="A") %>% filter(SampleExt %in% for_avg$Sample) %>% arrange(Lineage, phylo_order)

#average the reads
make_table <- data.frame(Region=ds$Region, Read=ds$Read, AAseq=ds$AAseq)

#i <- 5
for (i in 1:nrow(for_avg)) {
  o <- for_avg[i,1]
  oA <- paste0(o, "_A")
  oB <- paste0(o, "_B")
  df <- data.frame(A=ds[[oA]], B=ds[[oB]])
  df$avg <- (df$A + df$B) / 2
  df$avg_zeroed <- ifelse(df$A > error_threshold & df$B > error_threshold, df$avg, 0) 
  df <- df %>% select(avg_zeroed)
  colnames(df) <- o
  make_table <- cbind(make_table, df)
}

#filter to remove extra reads post zeroing:
keep2 <- rowSums(make_table[4:ncol(make_table)], na.rm=TRUE) > 0
sum(keep2)

final_data_seqs <- make_table[keep2,]


final_data_seqs <- final_data_seqs[,c("Region", "Read", "AAseq", sample_order$SampleExt)]
final_data_seqs <- final_data_seqs %>% arrange(Region, desc(across(all_of((sample_order$SampleExt)))))



write.csv(final_data_seqs, "~/Desktop/FINAL_avg_final_data_seqs.csv", row.names=FALSE)



#make percents
V1abs <- final_data_seqs %>% filter(Region=="V1")
V1perc <- V1abs %>% mutate_if(is.numeric, funs(./sum(.)*100)) 

V2abs <- final_data_seqs %>% filter(Region=="V2")
V2perc <- V2abs %>% mutate_if(is.numeric, funs(./sum(.)*100)) 

V3abs <- final_data_seqs %>% filter(Region=="V3")
V3perc <- V3abs %>% mutate_if(is.numeric, funs(./sum(.)*100)) 

V4abs <- final_data_seqs %>% filter(Region=="V4")
V4perc <- V4abs %>% mutate_if(is.numeric, funs(./sum(.)*100)) 

V5abs <- final_data_seqs %>% filter(Region=="V5")
V5perc <- V5abs %>% mutate_if(is.numeric, funs(./sum(.)*100)) 

V6abs <- final_data_seqs %>% filter(Region=="V6")
V6perc <- V6abs %>% mutate_if(is.numeric, funs(./sum(.)*100)) 

V7abs <- final_data_seqs %>% filter(Region=="V7")
V7perc <- V7abs %>% mutate_if(is.numeric, funs(./sum(.)*100)) 

all_perc <- rbind(V1perc, V2perc, V3perc, V4perc, V5perc, V6perc, V7perc)



#now make corr coeff plots to show similarity within clade and patient matched more similar: 

cor_plot <- all_perc[,4:ncol(all_perc)]
colnames(cor_plot) <- gsub("_extA", "", colnames(cor_plot))
colnames(cor_plot) <- gsub("_extB", "", colnames(cor_plot))
identical(colnames(cor_plot), sample_order$Sample) #TRUE, needs to be true. 

cor_df <- as.data.frame(matrix(nrow=ncol(cor_plot), ncol=ncol(cor_plot)))
colnames(cor_df) <- colnames(cor_plot)
rownames(cor_df) <- colnames(cor_plot)
pval_df <- cor_df

#i <- 5
for (i in 1:nrow(cor_df)) {
 # j <- 3
  for (j in 1:ncol(cor_plot)) {
  df2 <- data.frame(x=cor_plot[,i], y=cor_plot[,j])
  rm0 <- rowSums(df2) > 0 
  df2 <- df2[rm0,]
  hh <- cor.test(df2$x, df2$y)
  hh2 <- hh$estimate
  hh3 <- hh$p.value
  cor_df[i,j] <- hh2
  pval_df[i,j] <- hh3
  }
}

write.csv(cor_df, "~/Desktop/per_sample_tprk_pearsons.csv") 


