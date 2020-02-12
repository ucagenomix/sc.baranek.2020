library(ggplot2)
library(tidyr)


## Load / Create main variables and functions -----------------------------------------------------------

cell_type_order <- c('NKT0', 'Cycling NKT', 'NKT2a', 'NKT2b', 'NKT17', 'NKT1a', 'NKT1b', 'NKT1c', 'NKT1d')
method_color <-  c('wt' = '#009933', 'ko' = '#990000')

flexible_normalization <- function(data_in, by_row=TRUE){
  if(by_row){
    row_mean <- apply(data_in,1,mean)
    row_sd   <- apply(data_in,1,sd)
    output <- data_in
    for(i in 1:dim(data_in)[1]){
      output[i,] <- (data_in[i,] - row_mean[i])/row_sd[i]
    }
  }
  #### if by column
  if(!by_row){
    col_mean <- apply(data_in,2,mean)
    col_sd   <- apply(data_in,2,sd)
    output <- data_in
    for(i in 1:dim(data_in)[2]){
      output[,i] <- (data_in[,i] - col_mean[i])/col_sd[i]
    }
  }
  return(output)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

##### Heatmap WT ##########################################################################
setwd("/data/10x_data/000-notebooks/10x_paget/")
#zz=gzfile('/data/10x_data/000-notebooks/10x_paget/cellbrowser/wt/exprMatrix.tsv.gz', 'rt')
#countTable=read.table(zz, header=T, sep = "\t", row.names = 1, stringsAsFactors = F)
#save(countTable, file = '/data/10x_data/000-notebooks/10x_paget/cellbrowser/wt/exprMatrix.Rda')
load('/data/10x_data/000-notebooks/10x_paget/cellbrowser/wt/exprMatrix.Rda')

metadata <- read.table(file = "/data/10x_data/000-notebooks/10x_paget/cellbrowser/wt/meta.tsv",sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
coord <- read.table(file = "/data/10x_data/000-notebooks/10x_paget/cellbrowser/wt/umap_coords.tsv", sep = "\t", header = T)
cluster_color <- read.csv(file = "/data/10x_data/000-notebooks/10x_paget/cellbrowser/wt/clusterColor.tsv", sep = "\t", stringsAsFactors = F, header = T)
cell_type_color <- cluster_color$color
names(cell_type_color) <- cluster_color$Cluster
metadata$x = coord$x
metadata$y = coord$y
colnames(metadata)

rownames(countTable) <- sapply(strsplit(rownames(countTable), "[|]"), `[`, 2)

markers <- c('Lef1','Itm2a','Ccr9','Id3','Ldhb', # NKT0
             'Cks1b','Hmgb2','Mki67','Cenpa','Ran', # Cycling
             'Zbtb16','Plac8','Cmtm7','Icos','Izumo1r','Il1r1',  # NKT2
             'Rorc','Ccr6','Emb','Tmem176a', # NKT17
             'Ly6a','Fhl2','Itgae','Eng', # NKT1a
             'Id2','Nkg7','Klrb1c','Klrd1','Il2rb','Gimap6', # NKT1b
             'Klf2','S1pr1','Gzma','Fcer1g','Xcl1','S100a4', # NKT1c
             'Ifit3','Isg15','Stat1','Irf7','Bst2' # NKT1d
)

markers_flt_count <- countTable[markers, ]
colnames(markers_flt_count) <- gsub(".", "-", colnames(countTable), fixed = T)
colnames(countTable) <- gsub(".", "-", colnames(countTable), fixed = T)

markers_count <- matrix(0, nrow = length(markers), ncol = length(cell_type_order))
colnames(markers_count) <- cell_type_order
rownames(markers_count) <- markers

for (i in 1:ncol(markers_count)){
  print(colnames(markers_count)[i])
  cell_names <- rownames(metadata[metadata$cell_type == colnames(markers_count)[i], ])
  markers_count[, colnames(markers_count)[i]] <- rowMeans(markers_flt_count[, cell_names])
}


markers_count <- flexible_normalization(markers_count, by_row = T)
markers_count[markers_count > 2] <- 2
markers_count[markers_count < -1.5] <- -1.5

annotation <- data.frame(CellType = factor(colnames(markers_count), levels = cell_type_order))
rownames(annotation) <- colnames(markers_count)
ann_colors <- list(CellType = cell_type_color)

pp <- pheatmap::pheatmap(markers_count, cluster_rows = F, cluster_cols = F, 
                         annotation_col = annotation, annotation_colors = ann_colors,     
                         angle_col = 90, annotation_legend = F)

save_pheatmap_pdf(pp, "figures/Heatmap.wt.pdf", width = 7, height = 8)



##### # de ko-wt pop by pop ########## 
setwd("/data/10x_data/000-notebooks/10x_paget")
markers = read.delim("de.wt.ko.nkt17.csv", sep=",", stringsAsFactors = F, header=T)
head(markers)
#markers <- markers[abs(markers$pvals_adj) > 0.5,]

gg <- ggplot(data=markers, aes(x=logFC,y=pvals_adj)) + geom_point(color="darkgrey") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle(paste("NKT17 (ko-wt)", sep=""))

gg + geom_text_repel(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 3), aes(label=genes), size=3) +
     geom_point(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 3), col="red")

markers = read.delim("de.wt.ko.nkt2a.csv", sep=",", stringsAsFactors = F, header=T)
head(markers)
markers <- markers[abs(markers$logFC) < 4,]
gg <- ggplot(data=markers, aes(x=logFC,y=pvals_adj)) + geom_point(color="darkgrey") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle(paste("NKT2a (ko-wt)", sep=""))

gg + geom_text_repel(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 10), aes(label=genes), size=3) +
  geom_point(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 10), col="red")

markers = read.delim("de.wt.ko.nkt2b.csv", sep=",", stringsAsFactors = F, header=T)
head(markers)
markers <- markers[abs(markers$logFC) < 4,]
gg <- ggplot(data=markers, aes(x=logFC,y=pvals_adj)) + geom_point(color="darkgrey") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle(paste("NKT2b (ko-wt)", sep=""))

gg + geom_text_repel(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 10), aes(label=genes), size=3) +
  geom_point(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 10), col="red")

markers = read.delim("de.wt.ko.nkt1b.csv", sep=",", stringsAsFactors = F, header=T)
head(markers)
markers <- markers[abs(markers$logFC) < 4,]
gg <- ggplot(data=markers, aes(x=logFC,y=pvals_adj)) + geom_point(color="darkgrey") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle(paste("NKT1b (ko-wt)", sep=""))

gg + geom_text_repel(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 20), aes(label=genes), size=3) +
  geom_point(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 20), col="red")


markers = read.delim("de.wt.ko.nkt1c.csv", sep=",", stringsAsFactors = F, header=T)
head(markers)
markers <- markers[abs(markers$logFC) < 4,]
gg <- ggplot(data=markers, aes(x=logFC,y=pvals_adj)) + geom_point(color="darkgrey") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle(paste("NKT1c (ko-wt)", sep=""))

gg + geom_text_repel(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 10), aes(label=genes), size=3) +
  geom_point(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 10), col="red")

markers = read.delim("de.wt.ko.nkt1d.csv", sep=",", stringsAsFactors = F, header=T)
head(markers)
markers <- markers[abs(markers$logFC) < 4,]
gg <- ggplot(data=markers, aes(x=logFC,y=pvals_adj)) + geom_point(color="darkgrey") +
  theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank()) +
  ggtitle(paste("NKT1d (ko-wt)", sep=""))

gg + geom_text_repel(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 3), aes(label=genes), size=3) +
  geom_point(data=subset(markers, abs(markers$logFC) > 0.5 & markers$pvals_adj > 3), col="red")



