library(ggplot2)
library(tidyverse)
library(viridis)
library(ggdendro)
library(data.table)
library(scico)
library(gtable)
library(gridExtra)
library(cowplot)
library(RColorBrewer)
library(readxl)
library(xlsx)
library(ggsci)

# load DEXseq results
load("~/Desktop/sf3b1_splicing_figure/DEXSeq-AnalysisComplete.RData")

#### Create volcano plot ####

# create a df and drop rows with NA
dxr1_a <- as.data.frame(dxr1) %>% drop_na()
# filter by pvalue, log2 fold change, and add tags: up, down, or no sig change (but still keep unsig exons in df)
assign_sig <- function(df) {
  for (i in 1:nrow(df)) {
    print(i)
    if (df[i, "padj"] < 0.001) {
      if (df[i, "log2fold_sf3b1_control"] > 2) {
        df$significance[i] <- "Exon gained"
      } else if (df[i, "log2fold_sf3b1_control"] < -2) {
        df$significance[i] <- "Exon lost"
      }
    } else {
      df$significance[i] <- "No change"
    }
  }
  return(df)
}
# create df w/ filtering strategy
dxr1_df <- assign_sig(dxr1_a) 
dxr1_sig <- filter(dxr1_df, significance %in% c("Exon gained", "Exon lost"))
# remove Inf by adding back min padj value
min_padj <- min(filter(dxr1_df, padj != Inf & padj > 0)$padj)
# log transform padj
dxr1_df$padj1 <- log10(dxr1_df$padj + min_padj) * -1
# find min sig padj
min_sig <- min(filter(dxr1_df, significance %in% c("Exon gained", "Exon lost"))$padj1)
# create x axis ceiling
dxr1_df$log2fold_sf3b1_control[dxr1_df$log2fold_sf3b1_control > 20] <- 20
dxr1_df$log2fold_sf3b1_control[dxr1_df$log2fold_sf3b1_control < -20] <- -20
# plot
ggplot(data = dxr1_df, aes(x = log2fold_sf3b1_control, y = padj1)) + geom_bin2d(size =0.1, bins = 200) + scale_fill_continuous(type = "viridis", trans = "log2") + theme_bw() + ylab("-log10(padj)") + xlab("log2 fold change") + geom_vline(xintercept=c(-2,2), linetype = "dashed") + geom_hline(yintercept = min_sig, linetype = "dashed") + theme(text = element_text(size=16)) + labs(fill="Density bin count") + xlim(-20,20)


#### Create heatmap using normalized exon counts ####

# using counts in dxd
# create counts df
get_counts <- function(x, norm) {
  c <- as.data.frame(counts(x, normalized = norm))[,c(1:6)]
  names(c) <- c("control_rep1", "control_rep2", "control_rep3", "sf3b1_rep1", "sf3b1_rep2", "sf3b1_rep3")
  return(c)
}
counts <- get_counts(dxd,F)
counts_norm <- get_counts(dxd,T)
# add row mean and row sd and end of df rows
counts$Mean <- rowMeans(counts) 
counts.t <- transform(counts, Stdev=apply(counts[1:6], 1, sd))
# calculate COV (stdev/mean)
counts.t$COV <- counts.t$Stdev / counts.t$Mean
# indicate coverage filter by "pass" "fail"
counts.t$cov_filter <- apply(counts.t[1:6], 1, function(x) ifelse(test = length(which(x > 10)) > 2, yes = "pass", no = "fail"))
# filter by expression (at least 2 counts are > 10) and COV (> 0.7)
counts.t <- counts.t %>% filter(cov_filter == "pass" & COV > 0.7)
# subset normalized counts
pass_exons <- rownames(counts.t)
counts.n <- counts_norm %>% filter(rownames(counts_norm) %in% pass_exons)
# log2 transform norm counts
counts.nt <- log(counts.n + 1)
# format for clustering / plotting
counts.ng <-rownames_to_column(counts.nt, "Exons") %>% gather(key = "Condition", value = "Expression", -Exons) 
counts.s <- spread(counts.ng, key = Exons, value = Expression)
counts.s <- column_to_rownames(counts.s, "Condition")
# cluster
cluster_data <- function(x) {
  c <- x %>% dist(method="euclidean") %>% hclust(method = "ward.D2")
  dendro <- as.dendrogram(c)
  ddata <- dendro_data(dendro, type="rectangle")
  return(ddata)
}
exon_order <- cluster_data(counts.nt)
condition_order <- cluster_data(counts.s)
# order data
counts.ng$Exons <- factor(counts.ng$Exons, levels=as.character(exon_order[["labels"]]$label))
counts.ng$Condition <- factor(counts.ng$Condition, levels=as.character(condition_order[["labels"]]$label))
# color palette
cols <- brewer.pal(3, "YlOrRd")
pal <- colorRampPalette(cols)
# plot
heatmap <- ggplot(counts.ng, aes(x=Condition, y=Exons, fill=Expression)) + geom_tile() + scale_fill_gradientn(colors = pal(20), breaks=c(0,2,4,6,8)) + theme_bw() + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.text.x=element_text(angle=45, hjust=1), text = element_text(size=16)) + coord_cartesian(expand=F)

#### Create barplot from categorize_exons.py results ####

# read in df
exon_usage = read_csv("exon_usage_summary.csv", col_names = c("Exon Category", "Count"))
# modify names and order categories
exon_usage[exon_usage == 'KES'] <- 'Known ES'
exon_usage[exon_usage == 'NES'] <- 'Novel ES'
exon_usage$`Exon Category` <- factor(exon_usage$`Exon Category`, levels = c("TSS", "TTS", "Known ES", "Novel ES"))
# plot
ggplot(exon_usage, aes(`Exon Category`, Count, fill=`Exon Category`)) + geom_bar(stat='identity') + theme_bw() + theme(axis.text.x=element_text(angle=45, hjust=1), axis.title.x = element_blank(), legend.position = "none", text = element_text(size=16)) + coord_trans(y="sqrt") + scale_y_continuous(breaks=c(0,10,100,250,500,750,1000)) + scale_fill_jama()
