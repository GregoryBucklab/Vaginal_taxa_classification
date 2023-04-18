##### get human gene expression data (run two times to get pop2 and ptb47) ###########
setwd('/Users/BinZhu/secure/godel/gpfs_fs/bccl/bzhu/human_gene_expression/')
unique_gene = read.table('unique_gene.txt')
unique_gene = unique_gene$V1
unique_gene = unique_gene[order(unique_gene)]

setwd('/Users/BinZhu/secure/godel/gpfs_fs/bccl/bzhu/human_gene_expression/reads')
file_list = list.files(path = ".", pattern = '*.reads.txt')
sample = str_replace_all(file_list,'.reads.txt','')

reads_table_human_MTG = matrix(data = 0, ncol = length(file_list), nrow = length(unique_gene))
colnames(reads_table_human_MTG) = sample
row.names(reads_table_human_MTG) = unique_gene

for (a in 1 : length(file_list)) {
  setwd('/Users/BinZhu/secure/godel/gpfs_fs/bccl/bzhu/human_gene_expression/reads')
  reads_table = read.table(file_list[a])
  setwd('/Users/binzhu/Desktop/human_gene_expression/')
  gene_0 = unique_gene[!(unique_gene %in% reads_table$V1)]
  
  reads_table_2 = as.data.frame(matrix(data = NA, ncol = 2, nrow = length(gene_0)))
  colnames(reads_table_2) = c('V1','V2')
  reads_table_2$V1 = gene_0 
  reads_table_2$V2 = 0
  
  reads_table = rbind(reads_table,reads_table_2)
  reads_table=reads_table[order(reads_table$V1),]
  
  reads_table_human_MTG[,a] = reads_table$V2
}

reads_table_human_MTG = as.data.frame(reads_table_human_MTG)

setwd('/Users/binzhu/Desktop/human_gene_expression/')
reads_table_human_MTG_2 = as.data.frame(matrix(data = 0, ncol = ncol(reads_table_human_MTG)/2, nrow = nrow(reads_table_human_MTG)))
sample = str_replace_all(sample,'_.*','')
sample = sample[!duplicated(sample)]

colnames(reads_table_human_MTG) = sample
row.names(reads_table_human_MTG) = unique_gene

for (a in 1: length(sample)) {
  reads_table_human_MTG_2[,a] = reads_table_human_MTG[,(2*a-1)] + reads_table_human_MTG[,2*a]
}
reads_table_human_MTG = reads_table_human_MTG_2
rm(reads_table_human_MTG_2)

# gene reads threshold
keep = rowSums(reads_table_human_MTG) >= ncol(reads_table_human_MTG)
sum(keep)
reads_table_human_MTG = reads_table_human_MTG[keep,]

colnames(reads_table_human_MTG) = sample
row.names(reads_table_human_MTG) = unique_gene[keep]

colnames(reads_table_human_MTG) = str_replace_all(colnames(reads_table_human_MTG),'MV1R','MV1D')


# get annotation
gene_annotation = read.csv('human_gene_info.csv')

gene_annotation$order = gene_annotation$Type
gene_annotation$order[gene_annotation$order == 'mRNA'] = 'A'
gene_annotation$order[gene_annotation$order == 'snRNA'] = 'B'
gene_annotation$order[gene_annotation$order == 'lnc_RNA'] = 'C'
gene_annotation$order[gene_annotation$order == 'rRNA'] = 'D'
gene_annotation$order[gene_annotation$order == 'snoRNA'] = 'E'
gene_annotation$order[gene_annotation$order == 'guide_RNA'] = 'F'
gene_annotation$order[gene_annotation$order == 'scRNA'] = 'G'
gene_annotation$order[gene_annotation$order == 'antisense_RNA'] = 'H'
gene_annotation$order[gene_annotation$order == 'Y_RNA'] = 'I'
gene_annotation$order[gene_annotation$order == 'telomerase_RNA'] = 'J'
gene_annotation$order[gene_annotation$order == 'RNase_P_RNA'] = 'K'
gene_annotation$order[gene_annotation$order == 'RNase_MRP_RN'] = 'L'
gene_annotation$order[gene_annotation$order == 'transcript'] = 'M'
gene_annotation$order[gene_annotation$order == 'primary_transcript'] = 'N'

keep = gene_annotation$order %in% c('A','B','C','D','E','F','G','H','I','J','K','L','M','N')
gene_annotation = gene_annotation[keep,]

keep = row.names(reads_table_human_MTG) %in% gene_annotation$ID
sum(keep)
which(!keep)
reads_table_human_MTG = reads_table_human_MTG[-which(!keep),]

keep = gene_annotation$ID %in% row.names(reads_table_human_MTG)
sum(keep)
gene_annotation = gene_annotation[keep,]
gene_annotation = gene_annotation[order(gene_annotation$order),]
keep = !duplicated(gene_annotation$ID)
gene_annotation = gene_annotation[keep,]

#reads_table_human_MTG$name = gene_annotation$Gene

unique_gene_2 = unique(gene_annotation$Gene)
reads_table_human_MTG_2 = as.data.frame(matrix(data = NA, nrow = length(unique_gene_2), ncol = ncol(reads_table_human_MTG)))
colnames(reads_table_human_MTG_2) = colnames(reads_table_human_MTG)
row.names(reads_table_human_MTG_2) = unique_gene_2

for (a in 1: length(unique_gene_2)) {
  n = which(gene_annotation$Gene == unique_gene_2[a])
  reads_table_human_MTG_2[a,] = colSums(reads_table_human_MTG[n,])
}
reads_table_human_MTG = reads_table_human_MTG_2
rm(reads_table_human_MTG_2)

keep = !duplicated(gene_annotation$Gene)
gene_annotation = gene_annotation[keep,]



##### get 16s reads table #########
setwd('/Users/binzhu/secure/momspi-data/home/metadata/dataDrops/16s/')
data_1 = read.delim('pop2.stirrups.summary.20170606.txt')
data_2 = read.delim('ptb47.stirrups.summary.20170724.txt')

setwd('/Users/binzhu/Desktop/pH')

data = rbind(data_1,data_2)
data = data[data$SampleID %in% colnames(reads_table_human_MTG),]
data = data[data$Threshold.Status == 'AT',]

sample_list = unique(data$SampleID)

data$Taxa_2 = paste(data$Taxa, data$Threshold.Status, sep = '_')
data$Taxa_2 = str_replace_all(data$Taxa_2, '_AT','')
data$Taxa_2 = as.character(data$Taxa_2)

taxa_list = unique(data$Taxa_2)
taxa_list = taxa_list[order(taxa_list)]

trials = c(1: length(sample_list))

func_1 = function(trial) {
  keep = data$SampleID == sample_list[trial]
  data_2 = data[keep,]
  data_2 = data_2[,c(8,5)]
  keep = duplicated(data_2$Taxa_2)
  
  if (sum(keep) >0) {
    data_2 = data_2[order(data_2$Taxa_2),]
    keep = duplicated(data_2$Taxa_2)
    keep = which(keep)
    data_2$No_of_Reads[keep-1] = data_2$No_of_Reads[keep-1] + data_2$No_of_Reads[keep]
    data_2 = data_2[-keep,]
  }
  
  keep = taxa_list %in% data_2$Taxa_2
  sum(keep)
  data_3 = taxa_list[!keep]
  data_3 = as.data.frame(data_3)
  colnames(data_3) = 'Taxa_2'
  data_3$No_of_Reads = 0
  
  data_2 = rbind(data_2,data_3)
  data_2 = data_2[order(data_2$Taxa_2),]
  colnames(data_2) = c('Taxa',sample_list[trial])
  data_2$Taxa = NULL
  return(data_2)
}
data_3 = mclapply(trials, func_1, mc.cores = 8)

reads_table_16s = as.data.frame(matrix(data = 0, ncol = length(sample_list) , nrow = length(taxa_list)))
colnames(reads_table_16s) = sample_list
row.names(reads_table_16s) = taxa_list

for (a in 1: length(data_3)) {
  data_4 = data_3[a]
  data_3[a] = NA
  data_4 = unlist(data_4)
  reads_table_16s[,a] = as.numeric(as.character(data_4))
}
row.names(reads_table_16s)[str_detect(row.names(reads_table_16s), 'TM7_OTU-H1')] = "TM7_OTU_H1"

################## prepare reads tables ###################
{
  row.names(reads_table_16s)[row.names(reads_table_16s) == 'Lactobacillus_crispatus_cluster'] = 'Lactobacillus_crispatus'
  row.names(reads_table_16s)[row.names(reads_table_16s) == 'Lactobacillus_gasseri_cluster'] = 'Lactobacillus_gasseri'
  row.names(reads_table_16s)[row.names(reads_table_16s) == 'Lachnospiraceae_BVAB1'] = 'Ca._L._vaginae'
  row.names(reads_table_16s)[row.names(reads_table_16s) == 'TM7_OTU_H1'] = 'BVAB_TM7'
  
  taxa_list = read.csv('taxa_list.csv')
  taxa_list = taxa_list$x
  
  keep = str_detect(row.names(reads_table_16s), 'Lactobacillus|Gardnerella|Lachnospiraceae|Atopobium|Sneathia|Prevotella') |
    row.names(reads_table_16s) %in% taxa_list
  row.names(reads_table_16s)[keep]
  
  reads_table_1 = reads_table_16s[keep,]
  reads_table_2 = reads_table_16s[!keep,]

  keep = row.names(reads_table_1) %in% taxa_list
  sum(keep)
  reads_table_1_1 = reads_table_1[keep,]
  reads_table_1_2 = reads_table_1[!keep,]
  
  # reads_table_3
  taxa_list = str_replace_all(row.names(reads_table_2), '_.*', '')
  taxa_list = unique(taxa_list)
  
  reads_table_3 = as.data.frame(matrix(data = 0, ncol = ncol(reads_table_2), nrow = length(taxa_list)))
  row.names(reads_table_3) = taxa_list
  colnames(reads_table_3) = colnames(reads_table_2)
  
  taxa_list = str_replace_all(row.names(reads_table_2), '_.*', '')
  for (a in 1:nrow(reads_table_3)) {
    n = which(taxa_list == row.names(reads_table_3)[a])
    reads_table_3[a,] = colSums(reads_table_2[n,])
  }
  
  # reads_table_4
  taxa_list = str_replace_all(row.names(reads_table_1_2), '_.*', '')
  taxa_list = unique(taxa_list)
  
  reads_table_4 = as.data.frame(matrix(data = 0, ncol = ncol(reads_table_1_2), nrow = length(taxa_list)))
  row.names(reads_table_4) = taxa_list
  colnames(reads_table_4) = colnames(reads_table_1_2)
  
  taxa_list = str_replace_all(row.names(reads_table_1_2), '_.*', '')
  for (a in 1:nrow(reads_table_4)) {
    n = which(taxa_list == row.names(reads_table_4)[a])
    reads_table_4[a,] = colSums(reads_table_1_2[n,])
  }
  row.names(reads_table_4) = paste0('Other_',row.names(reads_table_4))
  
  reads_table = rbind(reads_table_1_1, reads_table_3,reads_table_4)
  
  
  taxa_list = read.csv('taxa_list.csv')
  taxa_list = taxa_list$x
  
  sum(row.names(reads_table) %in% taxa_list)
  taxa_list[!taxa_list %in% row.names(reads_table)]
  
  reads_table_16s = reads_table[row.names(reads_table) %in% taxa_list,]

}

keep = colSums(reads_table_16s) >= 5000
sum(keep)
reads_table_16s = reads_table_16s[,keep]

keep = colnames(reads_table_human_MTG) %in% colnames(reads_table_16s)
reads_table_human_MTG_16s = reads_table_human_MTG[,keep]
reads_table_16s = reads_table_16s[colnames(reads_table_human_MTG_16s)]
reads_table_16s_2 =reads_table_16s

metadata = metadata_all[metadata_all$SampleID %in% colnames(reads_table_16s),]
metadata_2 = metadata[match(colnames(reads_table_16s), metadata$SampleID),]
metadata_2$SampleID = paste0("Sample_",c(1:ncol(reads_table_16s)))

participant_list = unique(metadata_2$ParticipantID)
for (a in 1: length(participant_list)) {
  n = which(metadata_2$ParticipantID == participant_list[a])
  metadata_2$ParticipantID[n] = paste0("Participant_",a)
}

keep = !duplicated(metadata_2$ParticipantID)
sum(keep)
metadata_2 = metadata_2[keep,]
reads_table_human_MTG_16s = reads_table_human_MTG_16s[,keep]
reads_table_16s= reads_table_16s[,keep]

metadata_2$SampleID = paste0("Sample_",c(1:ncol(reads_table_16s)))
colnames(reads_table_human_MTG_16s) = paste0("Sample_",c(1:ncol(reads_table_human_MTG_16s)))
colnames(reads_table_16s) = paste0("Sample_",c(1:ncol(reads_table_16s)))



setwd('/Users/binzhu/Desktop/pH/Correlation_MTS')
write.csv(reads_table_16s,'reads_table_16s_MTS.csv', row.names = T)
write.csv(reads_table_human_MTG_16s,'reads_table_human_MTG_16s.csv', row.names = T)
write.csv(metadata_2,'metadata_2.csv', row.names = T)


##### normalize ############# 
# normalize data frames 
reads_table_16s = as.data.frame(t(reads_table_16s))
reads_table_16s = reads_table_16s + 0.5
reads_table_16s  <- clr(reads_table_16s)
reads_table_16s = as.data.frame(t(reads_table_16s))
reads_table_human_MTG_16s_2 = as.matrix(reads_table_human_MTG_16s)
library(DESeq2)
reads_table_human_MTG_16s_2 = varianceStabilizingTransformation(reads_table_human_MTG_16s_2) # sample in column    reads = FPKM*colSum/10e6
reads_table_human_MTG_16s_2 = as.data.frame(reads_table_human_MTG_16s_2)
rm(reads_table_cytokines,reads_table_16s_2,metadata_all,metadata_all_cha,metadata_all_num,reads_table_16s_abundance,reads_table_human_MTG,reads_table_human_MTG_16s)
rm(reads_table_16s_2, reads_table_human_MTG_16s,reads_table_human_MTG, reads_table,reads_table_1,reads_table_1_2,reads_table_2,reads_table_3,reads_table_4)
rm(reads_table_1_1, gene_annotation,data_3, taxonomy)
setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
save.image("~/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa/MTS_correlation.RData")
setwd('/Users/binzhu/Desktop/pH/Correlation_MTS')

##################### run in fenn ############
# $ cd /vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa/
# $ conda activate ALDEx2
# $ R

load("MTS_correlation.RData")
# correlation 
reads_table1 = reads_table_16s; reads_table2 = reads_table_human_MTG_16s_2
pvalue = 0.05;cor_parameter= 0;style = 1; bar_max = 2; bar_min = -2; pheatmap_fontsize = 5
treeheight = 50; alpha = 0.05; mc.cores =60

reads_table1 = as.matrix(reads_table1)
reads_table2 = as.matrix(reads_table2)

correlation_data_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = 0))
colnames(correlation_data_all) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
library(parallel)

try(for (a in 1: nrow(reads_table1)) {
  correlation_data = as.data.frame(matrix(data = NA, ncol = 5, nrow = nrow(reads_table2)))
  colnames(correlation_data) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
  correlation_data$Gene = row.names(reads_table1)[a]
  correlation_data$Taxa = row.names(reads_table2)
  
  trials = c(1: nrow(reads_table2))
  
  func_1 = function(trial) {
    data1 = as.numeric(reads_table1[a,])
    data2 = as.numeric(reads_table2[trial,])
    
    keep = (!is.na(data1)) & (!is.na(data2))
    data1 = data1[keep]
    data2 = data2[keep]
    
    if (length(data1) <8) {
      return(c(list(pvalue = NA), list(Rvalue = NA)))
    }
    
    data3 = cor.test(data1, data2, method = "spearman")
    pvalue = data3$p.value
    Rvalue = data3$estimate
    
    return(c(list(pvalue = pvalue), list(Rvalue = Rvalue)))
  }
  
  pRvalue = mclapply(trials, func_1, mc.cores = mc.cores)
  pRvalue = (unlist(pRvalue))
  keep = seq(1, length(pRvalue), 2)
  correlation_data$Pvalue = pRvalue[keep]
  
  keep = seq(2, length(pRvalue), 2)
  correlation_data$Rvalue = pRvalue[keep]
  
  if (a ==1) {
    correlation_data_all = correlation_data
  } else {
    correlation_data_all = rbind(correlation_data_all,correlation_data)
  }
})

data4 <- p.adjust(correlation_data_all$Pvalue, method = "BH")
#data4 <- adjust.p(correlation_data_all$Pvalue, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
#data4 <- data4$adjp
#correlation_data_all$adj_p <- data4$adjusted.p
correlation_data_all$adj_p <-data4

keep = correlation_data_all$adj_p <= pvalue
correlation_data_all = correlation_data_all[keep,]

keep = abs(correlation_data_all$Rvalue) >= cor_parameter
correlation_data_all_2 = correlation_data_all[keep,]

keep = !is.na(correlation_data_all_2$Gene)
correlation_data_all_2 = correlation_data_all_2[keep,]

# get correlation matrix
unique_gene = unique(correlation_data_all_2$Gene)
unique_taxa = unique(correlation_data_all_2$Taxa)

p.yes.rr = as.data.frame(matrix(data = 0, nrow = length(unique_gene), ncol = length(unique_taxa) ))
row.names(p.yes.rr) = unique_gene
colnames(p.yes.rr) = unique_taxa

for (a in 1: nrow(correlation_data_all_2)) {
  x = which(unique_gene == correlation_data_all_2$Gene[a])
  y = which(unique_taxa == correlation_data_all_2$Taxa[a])
  
  p.yes.rr[x,y]= correlation_data_all_2$Rvalue[a]
}

write.csv(p.yes.rr, 'p.yes.rr.csv', quote = F, row.names = T)

################# back to local computer ###########
setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu/vmb_taxa')
p.yes.rr = read.csv('p.yes.rr.csv', row.names = 1)

setwd('/Users/binzhu/Desktop/pH/Correlation_MTS')
p.yes.rr_cp = p.yes.rr


p.yes.rr= p.yes.rr_cp
library(igraph)

{
  paletteLength <- 50
  pvalue = 0.05;cor_parameter= 0;style = 1; pheatmap_fontsize = 5
  treeheight = 50; alpha = 0.05; mc.cores =60
  bar_max = max(p.yes.rr,na.rm = T)
  bar_min = min(p.yes.rr,na.rm = T)
  
  library(pheatmap)
  
  hclust_method = "ward.D"
  distance = "euclidean"
  cluster = pheatmap(p.yes.rr, method = hclust_method,clustering_distance_rows = distance,clustering_distance_cols = distance)
  
  taxa_cluster = sort(cutree(cluster$tree_row, k=10))
  taxa_cluster = as.data.frame(taxa_cluster)
  taxa_cluster$taxa = row.names(taxa_cluster)
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 7] = 6
  
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 3] = 11
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 5] = 12
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 8] = 13
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 4] = 14
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 2] = 15
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 1] = 16
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 6] = 17
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 9] = 18
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 10] = 19
  taxa_cluster$taxa_cluster = taxa_cluster$taxa_cluster -10
  
  cluster_color = cluster$tree_row$labels[cluster$tree_row$order]
  taxa_cluster = taxa_cluster[match(cluster_color,taxa_cluster$taxa),]
  write.csv(taxa_cluster,'Correlation_MTS_taxa_cluster.csv')  
  taxa_cluster$taxa = NULL
  taxa_cluster$taxa_cluster = as.factor(taxa_cluster$taxa_cluster)
  
  
  gene_cluster = sort(cutree(cluster$tree_col, k=7))
  gene_cluster = as.data.frame(gene_cluster)
  gene_cluster$gene = row.names(gene_cluster)
  cluster_color = cluster$tree_col$labels[cluster$tree_col$order]
  gene_cluster = gene_cluster[match(cluster_color,gene_cluster$gene),]
  gene_cluster$gene = factor(gene_cluster$gene, levels = gene_cluster$gene)
  gene_cluster$gene_cluster[gene_cluster$gene_cluster == 6] = 4
  gene_cluster$gene_cluster[gene_cluster$gene_cluster == 5] = 7
  
  gene_cluster$gene_cluster[gene_cluster$gene_cluster == 1] = 11
  gene_cluster$gene_cluster[gene_cluster$gene_cluster == 3] = 12
  gene_cluster$gene_cluster[gene_cluster$gene_cluster == 4] = 13
  gene_cluster$gene_cluster[gene_cluster$gene_cluster == 2] = 14
  gene_cluster$gene_cluster[gene_cluster$gene_cluster == 7] = 15
  gene_cluster$gene_cluster = gene_cluster$gene_cluster - 10
  write.csv(gene_cluster,'Correlation_MTS_gene_cluster.csv')
  gene_cluster$gene = NULL
  gene_cluster$gene_cluster = as.factor(gene_cluster$gene_cluster)
  
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
  
  library(RColorBrewer)
  cols <- colorRampPalette(brewer.pal(5, "Set1"))
  cols = cols(5)
  mycolors = list(taxa_cluster = c("1"="#FFDC00","2"="#49C800","3"="#04DDFF","4"="#aec7e8",
                                   "5"="#ff7f0e","6"="#d62728","7"="#FF00D1", "8"="#e377c2",
                                   "9"="#9C9146"),
                  gene_cluster = c("1"=cols[1],"2"=cols[2],"3"=cols[3],"4"=cols[4],
                                   "5"=cols[5]))
  
  cluster = pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                     treeheight_row = treeheight, treeheight_col = treeheight, method = hclust_method,
                     clustering_distance_rows = distance,
                     clustering_distance_cols = distance,show_colnames = F, 
                     annotation_col = gene_cluster, annotation_colors = mycolors,annotation_row = taxa_cluster)
  
  pdf("Correlation_MTS_5.pdf", width=8, height=5)
  print(cluster)
  dev.off()
  
  
  ### 
  
  taxa_cluster = sort(cutree(cluster$tree_row, k=10))
  taxa_cluster = as.data.frame(taxa_cluster)
  taxa_cluster$taxa = row.names(taxa_cluster)
  cluster_color = cluster$tree_row$labels[cluster$tree_row$order]
  taxa_cluster = taxa_cluster[match(cluster_color,taxa_cluster$taxa),]
  taxa_cluster$taxa = factor(taxa_cluster$taxa, levels = taxa_cluster$taxa)
  
  gene_cluster = sort(cutree(cluster$tree_col, k=150))
  gene_cluster = as.data.frame(gene_cluster)
  gene_cluster$gene = row.names(gene_cluster)
  cluster_color = cluster$tree_col$labels[cluster$tree_col$order]
  gene_cluster = gene_cluster[match(cluster_color,gene_cluster$gene),]
  gene_cluster$gene = factor(gene_cluster$gene, levels = gene_cluster$gene)
  
  gene_cluster$gene_cluster_new_order = NA
  gene_cluster_new_order = unique(gene_cluster$gene_cluster)
  for (a in 1: length(gene_cluster_new_order)) {
    n = which(gene_cluster$gene_cluster == gene_cluster_new_order[a])
    gene_cluster$gene_cluster_new_order[n] = a
  }
  
  p.yes.rr = p.yes.rr[as.character(gene_cluster$gene)]
  p.yes.rr = as.data.frame(t(p.yes.rr))
  p.yes.rr = p.yes.rr[as.character(taxa_cluster$taxa)]
  
  gene_cluster = cbind(gene_cluster, as.data.frame((p.yes.rr)))
  
  write.csv(gene_cluster,'Correlation_MTS_gene_cluster_150.csv')
  
  
  
  gene_cluster$gene %in% colnames(p.yes.rr)
  detach("package:igraph", unload = TRUE, force = T)
  
  # hclust_method = "ward.D"
  # distance = "euclidean"
  
  x = dist(p.yes.rr, method = distance, diag = FALSE, upper = T)
  x= as.matrix(x)
  write.csv(x, 'distance_in_host_interaction.csv')
}
