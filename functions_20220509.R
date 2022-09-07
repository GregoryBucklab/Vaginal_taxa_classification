######## library ########
library(parallel)
library(vegan) # for diversity
library(stringr) 
library(ggplot2)
library(Rtsne) # for t-SNE
library(ALDEx2) # for diffential abundance
# library(pracma)
# library(scales)
library(GUniFrac) # for 'Rarefy'
library(tidyr) # for 'gather'
# library(DiscriMiner)
library(compositions) # for 'clr'
library(reshape2) # for 'melt' in bar plot
library(cp4p) # for fdr

# for network
library(Hmisc)
library(Matrix)  
library(corrplot)
library(pheatmap)
library(RColorBrewer)
library(qgraph)


######### prepare reads table and metadata #########
prepare_reads_table = function(reads_table, metadata, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 1) {
  # total reads threshold
  keep <- colSums(reads_table) >= total_reads_threshold 
  reads_table <- reads_table[,keep]
  metadata <- metadata[keep,]
  
  # get abundance table
  trials = c(1: ncol(reads_table))
  func_1 = function(trial) {
    reads_table_abundance <- reads_table[,trial] / colSums(reads_table)[trial]
    return(reads_table_abundance)
  }
  reads_table_abundance = mclapply(trials, func_1, mc.cores = mc.cores)
  reads_table_abundance = unlist(reads_table_abundance)
  
  reads_table_abundance_2 = as.data.frame(matrix(data =0, ncol = ncol(reads_table), nrow = nrow(reads_table)))
  row.names(reads_table_abundance_2) <- row.names(reads_table)
  colnames(reads_table_abundance_2) <- colnames(reads_table)
  
  for (a in 1: ncol(reads_table)) {
    reads_table_abundance_2[,a] = reads_table_abundance[((a-1)*nrow(reads_table)+1):(a*nrow(reads_table))]
  }
  reads_table_abundance = reads_table_abundance_2
  
  # species threshold
  keep <- rep(T, nrow(reads_table_abundance))
  
  trials = c(1: nrow(reads_table_abundance))
  func_1 = function(trial) {
    c = sum(reads_table_abundance[trial,] >= species_threshold) / ncol(reads_table_abundance) >= 0.05        # input
    d = sum(reads_table_abundance[trial,] >= species_threshold/10) / ncol(reads_table_abundance) >= 0.15        # input
    keep = c|d
    return(keep)
  }
  keep = mclapply(trials, func_1, mc.cores = mc.cores)
  keep = unlist(keep)
  
  reads_table <- reads_table[keep,]
  
  # total reads threshold
  keep <- colSums(reads_table) >= total_reads_threshold 
  reads_table <- reads_table[,keep]
  metadata <- metadata[keep,]
  
  return(c(list(reads_table = reads_table),list(metadata = metadata)))
}





######### prepare reads table #########
prepare_reads_table_2 = function(reads_table, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 1) {
  # total reads threshold
  keep <- colSums(reads_table) >= total_reads_threshold 
  reads_table <- reads_table[,keep]
  
  # get abundance table
  trials = c(1: ncol(reads_table))
  func_1 = function(trial) {
    reads_table_abundance <- reads_table[,trial] / colSums(reads_table)[trial]
    return(reads_table_abundance)
  }
  reads_table_abundance = mclapply(trials, func_1, mc.cores = mc.cores)
  reads_table_abundance = unlist(reads_table_abundance)
  
  reads_table_abundance_2 = as.data.frame(matrix(data =0, ncol = ncol(reads_table), nrow = nrow(reads_table)))
  row.names(reads_table_abundance_2) <- row.names(reads_table)
  colnames(reads_table_abundance_2) <- colnames(reads_table)
  
  for (a in 1: ncol(reads_table)) {
    reads_table_abundance_2[,a] = reads_table_abundance[((a-1)*nrow(reads_table)+1):(a*nrow(reads_table))]
  }
  reads_table_abundance = reads_table_abundance_2
  
  # species threshold
  keep <- rep(T, nrow(reads_table_abundance))
  
  trials = c(1: nrow(reads_table_abundance))
  func_1 = function(trial) {
    c = sum(reads_table_abundance[trial,] >= species_threshold) / ncol(reads_table_abundance) >= 0.05        # input
    d = sum(reads_table_abundance[trial,] >= species_threshold/10) / ncol(reads_table_abundance) >= 0.15        # input
    keep = c|d
    return(keep)
  }
  keep = mclapply(trials, func_1, mc.cores = mc.cores)
  keep = unlist(keep)
  
  reads_table <- reads_table[keep,]
  
  # total reads threshold
  keep <- colSums(reads_table) >= total_reads_threshold 
  reads_table <- reads_table[,keep]
  
  return(c(list(reads_table = reads_table)))
}

######### vagitype #######################
vagitype <- function(reads_table, th = 0.3) {
 # reads_table=reads_table2
  
  reads_table_abundance = get_abundance_table(reads_table)
  
  vagitype_1 = matrix(data = NA, ncol =1, nrow = ncol(reads_table_abundance))
  colnames(vagitype_1) = 'Vagitype'
  
  for (a in 1:ncol(reads_table_abundance)) {
    if (max(reads_table_abundance[,a]) < th) {
      vagitype_1[a,1]= 'No type'
    } else {
      n= which(reads_table_abundance[,a] == max(reads_table_abundance[,a]), arr.ind=TRUE)
      vagitype_1[a,1]= row.names(reads_table_abundance)[n][1]
    }
  }
  
  return(vagitype_1)
}




######### get abundance table ################ input samples in columns ########
get_abundance_table <- function(reads_table) {
  reads_table_abundance = reads_table
  reads_table_abundance = sapply(1: ncol(reads_table), function(j) (reads_table_abundance[,j] <- reads_table[,j] / colSums(reads_table)[j] ))
  colnames(reads_table_abundance) = colnames(reads_table)
  row.names(reads_table_abundance) = row.names(reads_table)
  
  return(reads_table_abundance)
}

######### rarefaction ### reads_table_faction ### input samples in columns ######
reads_table_faction = function(reads_table, rarefy_to = NA) {
  reads_table <- t(reads_table)
  
  if (is.na(rarefy_to)) {
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  } else {
    reads_table = Rarefy(reads_table, depth = rarefy_to)
  }

  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(t(reads_table))
  return(reads_table)
}






######### alpha diversity ### alpha_diversity ### input samples in cols ###########
alpha_diversity = function(reads_table, metadata = NA, factor_name = NA, paired = F,order = NA, rarefy_to = NA) {
  
  if (is.na(metadata)[1]) {
    print('no metadata')
    return(NA)
  }
  
  if (is.na(factor_name)) {
    print('no factor name')
    return(NA)
  }
  
  # rarefy to normalize data
  reads_table = t(reads_table)
  
  if (is.na(rarefy_to)) {
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  } else {
    reads_table = Rarefy(reads_table, depth = rarefy_to)
  }
  
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  # calculate diversity
  alpha.shannon_diversity <- data.frame(diversity(reads_table))
  alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table))
  alpha.ovserved_OTU <- data.frame(colSums(t(reads_table) != 0))
  
  alpha = as.data.frame(matrix(data = NA,ncol=3,nrow = nrow(reads_table)))
  colnames(alpha) = c('alpha.shannon','alpha.evenness','alpha.ovserved_OTU')
  
  alpha$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table.
  alpha$alpha.evenness <- alpha.evenness$diversity.reads_table.
  alpha$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table.....0.
  
  metadata = cbind(metadata, alpha)
  
  colnames(metadata)[1] = 'factor'
  
  metadata = metadata[order(metadata$factor),]
  
  if (is.na(order)[1] ) {
    metadata$factor <- factor(metadata$factor , levels = unique(metadata$factor))
  } else {
    metadata$factor <- factor(metadata$factor , levels = order)
  }
  
  
  alpha.shannon = ggplot(metadata, aes(x=factor, y=alpha.shannon)) +
    geom_boxplot(aes(fill = factor),outlier.shape=NA) +
    labs(x = NULL, y = "Shannon index", fill=factor_name)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          axis.text.x=element_blank(),
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  
  alpha.evenness = ggplot(metadata, aes(x=factor, y=alpha.evenness)) +
    geom_boxplot(aes(fill = factor),outlier.shape=NA) +
    labs(x = NULL, y = "Evenness", fill=factor_name)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          axis.text.x=element_blank(),
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  
  alpha.ovserved_OTU = ggplot(metadata, aes(x=factor, y=alpha.ovserved_OTU)) +
    geom_boxplot(aes(fill = factor),outlier.shape=NA) +
    labs(x = NULL, y = "Observed OTU", fill=factor_name)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          axis.text.x=element_blank(),
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  # geom_jitter(shape=16, position=position_jitter(0.2), size = 0.5)+
  
  # calculate significance
  factor_levels = unique(metadata$factor)
  n = length(factor_levels)
  
  Shannon_sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
  colnames(Shannon_sig) = factor_levels
  row.names(Shannon_sig) = factor_levels
  
  Evenness_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(Evenness_sig) = factor_levels
  row.names(Evenness_sig) = factor_levels
  
  OTU_sig = as.data.frame(matrix(data = NA, ncol=n, nrow = n))
  colnames(OTU_sig) = factor_levels
  row.names(OTU_sig) = factor_levels
  
  for (a in 1:(n-1)) {
    for (b in (a+1) : n) {
      factor_level1 <- subset(metadata,  factor == factor_levels[a],
                              drop = TRUE)
      factor_level2 <- subset(metadata,  factor == factor_levels[b],
                              drop = TRUE)
      
      Shannon_sig[a,b] <- wilcox.test(factor_level1$alpha.shannon, 
                                      factor_level2$alpha.shannon, paired = paired)$p.value
      Evenness_sig[a,b] <- wilcox.test(factor_level1$alpha.evenness, 
                                       factor_level2$alpha.evenness, paired = paired)$p.value
      OTU_sig[a,b] <- wilcox.test(factor_level1$alpha.ovserved_OTU, 
                                  factor_level2$alpha.ovserved_OTU, paired = paired)$p.value
      
    }
  }
  output = c(list(alpha = alpha), list(shannon = alpha.shannon), 
             list(evenness =alpha.evenness) , list(ovserved_OTU =alpha.ovserved_OTU),
             list(sig_Shannon = Shannon_sig),list(sig_Evenness = Evenness_sig),
             list(sig_OTU = OTU_sig))
  
  return(output)
}


######### beta diversity ### beta_diversity ### input samples in cols; metadata and factor_name are needed; order of factors could be set; ref_group is for setting the reference of bc distance in different groups; can skip from NMDS; output bc distance, within sample distance, distance amoung groups and NMDS #####
beta_diversity = function(reads_table, metadata = NA, factor_name = NA, order = NA, NMDS_skip = T, ref_group = NA, rarefy_to = NA, pheatmap_fontsize = 50,treeheight = 10, pheatmap_y = T) {
  
  if (is.na(metadata)[1]) {
    print('no metadata')
    return(NA)
  }
  
  if (is.na(factor_name)[1]) {
    print('no factor name')
    return(NA)
  }
  
  # rarefy to normalize data
  reads_table = as.data.frame(t(reads_table))
  
  if (is.na(rarefy_to)) {
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  } else {
    reads_table = Rarefy(reads_table, depth = rarefy_to)
  }
  
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  metadata=as.matrix(metadata)
  
  # Bray_Curtis
  Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
  Bray_Curtis <- as.data.frame(Bray_Curtis)
  
  Bray_Curtis_2 = Bray_Curtis
  #    Bray_Curtis_2[row(Bray_Curtis_2) <= col(Bray_Curtis_2)] =NA
  
  # within sample distance
  group_dis = gather(Bray_Curtis_2)
  group_dis$key2 = rep(row.names(Bray_Curtis_2),ncol(Bray_Curtis_2))
  
  Source = matrix(data = NA, ncol = length(metadata), nrow = length(metadata))
  
  for (a in 1:length(metadata)) {
    Source[a,] = metadata
  }
  Source = gather(as.data.frame(Source))
  group_dis$Source = Source$value
  group_dis$Target = rep(metadata,length(metadata))
  
  group_dis = group_dis[!is.na(group_dis$value),]
  group_dis$Source <- as.factor(group_dis$Source)
  #  group_dis$value = as.numeric(as.character(group_dis$value))
  
  keep = group_dis$Source == group_dis$Target
  within_dis = group_dis[keep,]
  keep = within_dis$key != within_dis$key2
  within_dis = within_dis[keep,]
  #  within_dis$value = as.numeric(as.character(within_dis$value))
  
  if (!is.na(order)) {
    within_dis$Source = as.factor(within_dis$Source)
    within_dis$Source = factor(within_dis$Source, levels= order)
  }
  
  within_dis_p = ggplot(within_dis, aes(x=Source, y=value)) +
    geom_boxplot(aes(fill = Source),outlier.shape=NA) +
    labs(x = NULL, y = "Within sample distance", fill=factor_name)+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
  within_dis_p
  ggsave("within_group_distance.pdf", width=3, height=2)
  
  
  # significance among within sample distance, Wilcoxon test
  group_level = unique(within_dis$Source)
  n = length(group_level)
  
  within_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
  colnames(within_dis_sig) = group_level
  row.names(within_dis_sig) = group_level
  for (a in 1:(n-1)) {
    for (b in (a+1): n) {
      set1 <- as.matrix(subset(within_dis,  Source == group_level[a], value,
                               drop = F))
      set2 <- as.matrix(subset(within_dis,  Source == group_level[b], value,
                               drop = F))
      within_dis_sig[a,b] <- wilcox.test(set1, set2, paired = F)$p.value
      
    }
  }
  write.csv(within_dis_sig,'within_group_distance.csv', row.names = T, quote = F)
  
  # distance among groups
  {
    group_level = unique(group_dis$Source)
    n = length(group_level)
    distance_median = matrix(data=NA, nrow = n, ncol =n)
    colnames(distance_median) = group_level
    row.names(distance_median) = group_level
    
    group_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
    colnames(group_dis_sig) = group_level
    row.names(group_dis_sig) = group_level
    
    for (a in 1:n) {
      for (b in a:n) {
        distance_data = group_dis$value[group_dis$Source == row.names(distance_median)[a] & group_dis$Target == colnames(distance_median)[b]]
        distance_median[a,b] <- median(distance_data)
        distance_median[b,a] <- distance_median[a,b]
        
        if (a != b) {
          keep = metadata == group_level[a] | metadata == group_level[b]
          metadata_2 = as.character(metadata[keep])
          reads_table_2 = reads_table[keep,]
          
          metadata_2 = as.data.frame(metadata_2)
          pvalue <- adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")[1,5] 
          group_dis_sig[b,a] = pvalue
          group_dis_sig[a,b] = pvalue
        }
        
      }
    }
    write.csv(group_dis_sig,'between_group_distance.csv', row.names = T, quote = F)
    
    if (pheatmap_y == F & n > 2) {
      group_dis_2_p = corrplot(distance_median,p.mat = group_dis_sig, method = 'shade', diag = F, type="upper",
                               sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
                               is.corr=FALSE, col=colorRampPalette(c("white","red"))(100)[c(10:100)],
                               order="original",tl.col = "black")
      
      pdf("between_group_distance_2.pdf")
      corrplot(distance_median,p.mat = group_dis_sig, method = 'shade', diag = F, type="upper",
               sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
               is.corr=FALSE, col=colorRampPalette(c("white","red"))(100)[c(10:100)],
               order="original",tl.col = "black")
      dev.off()
      
    } else {
      group_dis_sig_2 = group_dis_sig
      group_dis_sig_2[is.na(group_dis_sig)] = ''
      group_dis_sig_2[group_dis_sig > 0.05] = ''
      group_dis_sig_2[group_dis_sig <= 0.05 & group_dis_sig > 0.01] = '*'
      group_dis_sig_2[group_dis_sig <= 0.01 & group_dis_sig > 0.001] = '**'
      group_dis_sig_2[group_dis_sig <= 0.001] = '***'
      
      save_heatmap_pdf <- function(x, filename, width=8, height=8) {
        stopifnot(!missing(x))
        stopifnot(!missing(filename))
        pdf(filename, width=width, height=height)
        grid::grid.newpage()
        grid::grid.draw(x$gtable)
        dev.off()
      }
      
      group_dis_2_p = pheatmap(distance_median, cluster_rows=TRUE, show_rownames=TRUE, 
                               cluster_cols=T, show_colnames=T, 
                               color=colorRampPalette(c("white","red"))(100),
                               fontsize = pheatmap_fontsize, display_numbers = group_dis_sig_2, 
                               treeheight_row = treeheight, treeheight_col = treeheight)
      save_heatmap_pdf(group_dis_2_p, "between_group_distance_2.pdf", width=2, height=2)
    }
    
  }
  
  
  
  # Running Nonmetric Multidimensional Scaling (NMDS) Ordination
  if (NMDS_skip == T) {
    # output
    output = c(Bray_Curtis = list(Bray_Curtis), within_dis_p = list(within_dis_p),
               within_dis_sig = list(within_dis_sig), group_dis_2_p = list(group_dis_2_p),
               group_dis_sig = list(group_dis_sig))
    return(output)
    
  } else {
    
    colnames(metadata)[1] = 'factor'
    
    NMDS <-
      metaMDS(Bray_Curtis,
              distance = "bray",
              k = 2,
              maxit = 999, 
              trymax = 20,
              wascores = TRUE)
    
    mds_data <- as.data.frame(NMDS$points)
    mds_data$factor <- metadata
    
    if (!is.na(order)[1]) {
      mds_data$factor = as.factor(mds_data$factor)
      mds_data$factor = factor(mds_data$factor, levels= order)
    }
    
    NMDS = ggplot(mds_data, aes(x = MDS1, y = MDS2, color = factor)) +
      geom_point(size = 0.3)+
      scale_colour_discrete(factor_name)+
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    NMDS
    ggsave("NMDS.pdf", width=3, height=3)
    
    NMDS_2 = ggplot(mds_data, aes(x = MDS1, y = MDS2, color = factor)) +
      geom_point(size = 0.3)+
      scale_colour_discrete(factor_name)+
      stat_ellipse(type = "t")+
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    NMDS_2
    ggsave("NMDS_2.pdf", width=3, height=3)
    
    # output
    output = c(Bray_Curtis = list(Bray_Curtis), within_dis_p = list(within_dis_p), 
               within_dis_sig = list(within_dis_sig),group_dis_2_p = list(group_dis_2_p),
               group_dis_sig = list(group_dis_sig), NMDS =list(NMDS), NMDS_2 =list(NMDS_2))
    return(output)
  }
  
}





######### diffential abundance ### dif_abundance; sample in column; paired text or not; change order of species; fold change threshold; if style = 1, give dot plot, else give box plot ####################
dif_abundance <- function(reads_table,metadata, pvalue_th = 0.05, fold_change_th = 1, paired_test = F, order_reverse = F, style = 1, order = NA) {
  # style =1
  # reads_table= cts
  # metadata=as.character(metadata_virus_PTB$Delivery)
  # paired_test = F
  # order_reverse = F
  # style =2 
  # order = c('TB','PTB')
  # fold_change_th = 1
  
  conds <- metadata
  
  x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
  
  # paired Wilcoxon Rank Sum test and Welch's t-test
  x.tt <- aldex.ttest(x, paired.test= paired_test)
  
  if (min(x.tt$wi.eBH) > 0.05) {
    print('No taxon has significant abundance change')
    return(NA)
  }
  
  x.effect <- aldex.effect(x)
  
  x.all <- data.frame(cbind(x.tt,x.effect))
  
  abundance_change_das <- x.all$diff.btw # diff.btw is a vector containing the per-feature median difference between condition A and B
  
  if (max(abs(abundance_change_das)) < fold_change_th) {
    print('No taxon has significant abundance change')
    return(NA)
  }
  
  if (order_reverse == T) {
    abundance_change_das = -abundance_change_das
  }
  
  x.all <- cbind(abundance_change_das,x.all)
  
  `-Log10(adj-pvalue)` <- -log10(x.all$wi.eBH)     # use wi.eBH as the adj-pvalue
  #x.all$abundance_das <- abundance_das
  x.all$`-Log10(adj-pvalue)` <- `-Log10(adj-pvalue)`
  
  # draw figure
  das <- x.all[(x.all$`-Log10(adj-pvalue)` >= -log10(pvalue_th) & (x.all$abundance_change_das >=fold_change_th | x.all$abundance_change_das <=-fold_change_th)),]
  
  if (nrow(das)==0) {
    print('No taxon has significant abundance change')
    return(NA)
  }
  
  das$Species <- row.names(das)
  das <- das[order(das$abundance_change_das),] 
  
  metadata = as.factor(metadata)
  lev = levels(metadata)
  
  if (order_reverse == T) {
    das$Color <- ifelse(das$abundance_change_das < 0, paste0("Enriched in ",lev[2]), paste0("Enriched in ",lev[1]))  # above / below avg flag
  } else {
    das$Color <- ifelse(das$abundance_change_das < 0, paste0("Enriched in ",lev[1]), paste0("Enriched in ",lev[2]))  # above / below avg flag
    
  }
  
  das$Species <- factor(das$Species, levels = das$Species)  # convert to factor to retain sorted order in plot.
  
  das$abundance_change_das[das$abundance_change_das == Inf] = 10
  das$abundance_change_das[das$abundance_change_das == -Inf] = -10
  das$abundance_change_das[das$abundance_change_das <= -10] = -10
  das$abundance_change_das[das$abundance_change_das >= 10] = 10
  
  
  if (style == 1) {
    theme_set(theme_bw())  
    
    p <- ggplot(das, aes(Species, abundance_change_das)) + 
      geom_point(aes(col=Color, size=`-Log10(adj-pvalue)`)) + 
      coord_flip() +          # convert x y axis
      labs(x = 'Taxa', y = "Median difference in clr values")+ 
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))  
    
    
  } else {
    taxa_list = row.names(das)
    
    keep = which(row.names(reads_table) %in% taxa_list)
    
    reads_table_abundance = get_abundance_table(reads_table, mc.cores = 8)
    reads_table_2 <- reads_table_abundance[keep,]
    reads_table_2 <- as.data.frame(t(reads_table_2))
    
    #   colnames(reads_table_2)=taxa_list
    reads_table_3 = gather(reads_table_2)
    reads_table_3$Type = rep(conds, length(taxa_list))
    colnames(reads_table_3) = c('Taxa', 'Abundance','Type')
    
    reads_table_3$Taxa = str_replace_all(reads_table_3$Taxa,'.*__','')
    
    if (!is.na(order)[1]) {
      reads_table_3$Type = factor(reads_table_3$Type, levels = order)
    }
    
    p <- ggplot(reads_table_3, aes(x=Taxa, y=Abundance,fill=Type)) +
      geom_boxplot(outlier.shape = NA)+
      ylab("Abundance (%)")+
      theme(axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))
  }
  return(c(list(p = p), list(data = x.all)))
  
}


######### bar plot ########## samples on columns ###########
barplot <- function(reads_table, type_th = 0.3, taxa_num = 9) { # pass proportion data, samples on columns, taxa on rows.
  #   reads_table = reads_table_virus_2
  #   type_th = 0.1
  
  reads_table_abundance <- matrix(data =0, ncol = ncol(reads_table),nrow = nrow(reads_table))
  
  for (a in 1:ncol(reads_table)) {
    reads_table_abundance[,a] <- reads_table[,a] / colSums(reads_table)[a]
    
  }
  row.names(reads_table_abundance) = row.names(reads_table)
  colnames(reads_table_abundance) = colnames(reads_table) 
  
  reads_table = as.data.frame(t(reads_table_abundance))
  
  reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
  
  # assign vagitype
  mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
  maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
  mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
  mytypes[maxprop < type_th] <- "No Type"
  uniqueTypes <- unique(mytypes)
  
  top_species = as.data.frame(sort(colSums(reads_table),decreasing=T))
  
  keep = mytypes %in% row.names(top_species)[1:taxa_num]
  
  mytypes[!keep] = 'Others'
  
  if (length(uniqueTypes) < (taxa_num+1)) {
    myTypeOrder <- c( row.names(top_species)[1:(length(uniqueTypes)-1)], "No Type")
  } else {
    myTypeOrder <- c( row.names(top_species)[1:taxa_num], "No Type","Others")
  }
  myTypeOrder <- data.frame(typeOrder=c(1:length(myTypeOrder)),row.names=myTypeOrder)
  
  reads_table <- 100*reads_table
  
  # order data by type then proportion
  reads_table <- reads_table[order(myTypeOrder[mytypes,],-maxprop),]  # reorder the samples by vagitype and the abundance of dominant taxa
  barplottypes <- mytypes[order(myTypeOrder[mytypes,],-maxprop)]   # get the list of vagitypes of samples
  
  # color
  mycolors <- rep("gray", ncol(reads_table))
  
  color_code_list = c("#3264B8","#FC9800","#C80B0B","#FDFFBA","#55ff42",
                      "#F7A0A0","#D7A7F1","#337F08","#00FFE8","#9C58FE","#A68300",'gray')
  
  color_code_list = color_code_list[1:(taxa_num+1)]
  
  for (a in 1:(length(color_code_list)-1)) {
    mycolors[which(colnames(reads_table) == row.names(top_species)[a])] <- color_code_list[a]
  }
  
  mycolors <- factor(mycolors)
  mycolordf <- data.frame(mycolors, Taxa=colnames(reads_table))
  
  mycolordf2 = mycolordf[mycolordf$mycolors != 'gray',]
  mycolordf2 = mycolordf2[1:taxa_num,]
  mycolordf2[(taxa_num+1),1] = "gray"
  mycolordf2[(taxa_num+1),2] = "Others"
  mycolors_label = mycolordf2$Taxa[order(mycolordf2$mycolors)]
  
  sampleorder <- factor(rownames(reads_table))
  
  ggplotdata <- melt(as.matrix(reads_table),  
                     id.vars=names(reads_table), 
                     varnames=c("SampleID", "Taxa"), 
                     value.name="ATprop")
  ggplotdata$SampleID <- factor(ggplotdata$SampleID, levels=sampleorder)
  ggplotdata <- merge(ggplotdata, mycolordf, by="Taxa")
  ggplotdata <- ggplotdata[ggplotdata$ATprop != 0,]
  
  p <- ggplot(ggplotdata, aes(SampleID, ATprop, fill=mycolors, group=ATprop)) + 
    geom_bar(stat="identity", position="stack", width=1) +
    scale_fill_manual(values = levels(mycolors),labels = mycolors_label) +
    labs(x="Sample", y="Relative Abundance") + 
    theme(axis.text.x=element_blank(),
          axis.ticks=element_blank())  
  #,labels = mycolors_label
  output = c(p = list(p), color = list(mycolordf))
  return(output)
}




######### network rcorr ############# sample in cols #######
newwork_rcorr <- function(reads_table, pvalue = 0.05, cor_parameter= 0, style = 1, bar_max = 2, bar_min = -2, pheatmap_fontsize = 5, treeheight = 50, alpha = 0.05) {
  #  reads_table = reads_table_2
  
  reads_table = as.data.frame(t(reads_table))
  reads_table = reads_table + 0.5
  reads_table <- clr(reads_table)      ### CLR normalization in rows
  
  # pvalue is the pvalue threshold for cooccurence
  # cor_parameter is the cooccurence value threshold
  # style is the style of the heatmap
  
  #reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  #reads_table <- reads_table$otu.tab.rff
  #reads_table <- as.data.frame(reads_table)
  
  reads_table = as.matrix((reads_table))
  
  otu.cor <- rcorr(reads_table, type="spearman")
  
  otu.pval <- forceSymmetric(otu.cor$P)
  otu.pval <- otu.pval@x
  otu.pval <- adjust.p(otu.pval, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  otu.pval <- otu.pval$adjp
  otu.pval <- otu.pval$adjusted.p
  
  p.yes <- otu.pval< pvalue  
  
  r.val = otu.cor$r # select all the correlation values 
  p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion 
  
  p.yes.r <- abs(p.yes.r) > cor_parameter # output is logical vector
  p.yes.rr <- p.yes.r*r.val # use logical vector for subscripting.
  
  p.yes.rr[is.na(p.yes.rr)] = 0
  
  keep = abs(colSums(p.yes.rr)) > 0
  p.yes.rr = p.yes.rr[,keep]
  p.yes.rr = p.yes.rr[keep,]
  
  # gephi output
  gephi_p.yes.rr = as.matrix(p.yes.rr)
  for (a in 1:nrow(p.yes.rr)) {
    for (b in 1:ncol(p.yes.rr)) {
      if (a >=b) {
        gephi_p.yes.rr[a,b] = NA
      }
    }
  }
  gephi_p.yes.rr <- as.data.frame(gephi_p.yes.rr)
  gephi_p.yes.rr = gather(gephi_p.yes.rr)
  gephi_p.yes.rr$Taxa = rep(row.names(p.yes.rr),ncol(p.yes.rr))
  gephi_p.yes.rr = gephi_p.yes.rr[!is.na(gephi_p.yes.rr[,2]),]
  gephi_p.yes.rr$value = as.numeric(as.character(gephi_p.yes.rr$value))
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr$value != 0,] 
  gephi_p.yes.rr = gephi_p.yes.rr[gephi_p.yes.rr[,2] >= cor_parameter | gephi_p.yes.rr[,2] <= -cor_parameter,]
  colnames(gephi_p.yes.rr) = c('Source','Weight','Target')
  # heatmap
  library(igraph)
  
  if (style == 1) {
    if (bar_max == 2) {
      bar_max = max(p.yes.rr,na.rm = T)
      bar_min = min(p.yes.rr,na.rm = T)
      paletteLength <- 50
      myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
      myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                    seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
      
      myBreaks <- unique(myBreaks)
      
      if (sum(myBreaks < 0) == 0) {
        myColor <- colorRampPalette(c("white", "red"))(paletteLength)
        p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      } else if (sum(myBreaks > 0) == 0) {
        myColor <- colorRampPalette(c("blue", "white"))(paletteLength)
        p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      } else {
        p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      }
    } else {
      paletteLength <- 50
      myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
      myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                    seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
      
      if (sum(myBreaks < 0) == 0) {
        myColor <- colorRampPalette(c("white", "red"))(paletteLength)
        p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      } else if (sum(myBreaks > 0) == 0) {
        myColor <- colorRampPalette(c("blue", "white"))(paletteLength)
        p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      } else {
        p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      }
    }
    
  } else {
    p.yes.rr <- as.matrix(p.yes.rr)
    p = corrplot(p.yes.rr, type="upper", order="hclust", insig = "blank",tl.col = "black",na.label = "o")
  }
  
  detach("package:igraph", unload = TRUE, force = T)
  return(c(list(p = p), list(gephi_input = gephi_p.yes.rr), list(cor_matrix = p.yes.rr)))
}





######### newwork_rcorr3 correlation of two paired normalized reads table ############# sample in col #######
newwork_rcorr3 <- function(reads_table1, reads_table2, pvalue = 0.05, cor_parameter= 0, style = 1, bar_max = 2, bar_min = -2, pheatmap_fontsize = 5 , treeheight = 50, alpha = 0.05, mc.cores =1) {
  #    cor_parameter= 0
  #    pvalue = 0.05
  #    style = 1
  #    bar_max = 2
  #    bar_min = -2
  #    pheatmap_fontsize = 7 
  #    alpha = 0.05
  #    reads_table1=reads_table_human_MTG_16s
  #    reads_table2=alpha.shannon_diversity
  #    mc.cores = 8
  
  reads_table1 = as.matrix(reads_table1)
  reads_table2 = as.matrix(reads_table2)
  
  correlation_data_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = 0))
  colnames(correlation_data_all) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
  
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
      
      data3 = cor.test(data1, data2, method = "pearson")
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
  
  data4 <- adjust.p(correlation_data_all$Pvalue, pi0.method="bky", alpha = alpha,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  data4 <- data4$adjp
  correlation_data_all$adj_p <- data4$adjusted.p
  
  
  keep = correlation_data_all$adj_p <= pvalue
  correlation_data_all = correlation_data_all[keep,]
  
  keep = abs(correlation_data_all$Rvalue) >= cor_parameter
  correlation_data_all_2 = correlation_data_all[keep,]
  
  keep = !is.na(correlation_data_all_2$Gene)
  correlation_data_all_2 = correlation_data_all_2[keep,]
  
  if (nrow(correlation_data_all_2) == 0) {
    return(NA)
  }
  
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
  
  # heatmap
  if (length(p.yes.rr) == 0) {
    return(NA)
  }
  
  if (nrow(p.yes.rr) == 0 | ncol(p.yes.rr) == 0) {
    return(NA)
  }
  
  if (nrow(p.yes.rr) == 1 | ncol(p.yes.rr) == 1) {
    return(p.yes.rr)
  }
  
  library(igraph)
  
  if (style == 1) {
    if (bar_max == 2) {
      bar_max = max(p.yes.rr,na.rm = T)
      bar_min = min(p.yes.rr,na.rm = T)
      paletteLength <- 50
      
      if (bar_max > 0 & bar_min < 0) {
        myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
        myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                      seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
        p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
        
      } else if (bar_max > 0 & bar_min == 0) {
        myColor <- colorRampPalette(c("white", "red"))(paletteLength)
        p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      } else {
        myColor <- colorRampPalette(c("blue","white"))(paletteLength)
        p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      }
      
    } else {
      paletteLength <- 50
      myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
      myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                    seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
      
      p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                    treeheight_row = treeheight, treeheight_col = treeheight)
    }
    
  } else {
    p.yes.rr <- as.matrix(p.yes.rr)
    p = corrplot(p.yes.rr, type="upper", order="hclust", insig = "blank",tl.col = "black",na.label = "o")
  }
  
  detach("package:igraph", unload = TRUE, force = T)
  return(c(list(p = p), list(cor_matrix = p.yes.rr), list(cor_matrix2 = correlation_data_all_2)))
}



######### newwork_rcorr2 correlation of two paired normalized reads table - no FDR ############# sample in col #######
newwork_rcorr2 <- function(reads_table1, reads_table2, pvalue = 0.05, cor_parameter= 0, style = 1, bar_max = 2, bar_min = -2, pheatmap_fontsize = 5, treeheight = 50, alpha = 0.05, mc.cores =1) {
  #      cor_parameter= 0
  #      pvalue = 0.05
  #      style = 1
  #      bar_max = 2
  #      bar_min = -2
  #      pheatmap_fontsize = 7 
  #      alpha = 0.05
  #      reads_table1=reads_table2
  #      reads_table2=reads_table_4
  #      mc.cores = 8
  
  reads_table1 = as.matrix(reads_table1)
  reads_table2 = as.matrix(reads_table2)
  
  correlation_data_all = as.data.frame(matrix(data = NA, ncol = 5, nrow = 0))
  colnames(correlation_data_all) = c('Gene','Taxa','Pvalue','Rvalue','adj_p')
  
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
      
      data3 = cor.test(data1, data2, method = "pearson")
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
  
  keep = correlation_data_all$Pvalue <= pvalue
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
  
  # heatmap
  if (length(p.yes.rr) == 0) {
    return(NA)
  }
  
  if (nrow(p.yes.rr) == 0 | ncol(p.yes.rr) == 0) {
    return(NA)
  }
  
  if (nrow(p.yes.rr) == 1 | ncol(p.yes.rr) == 1) {
    return(p.yes.rr)
  }
  
  library(igraph)
  if (style == 1) {
    if (bar_max == 2) {
      bar_max = max(p.yes.rr,na.rm = T)
      bar_min = min(p.yes.rr,na.rm = T)
      paletteLength <- 50
      
      if (bar_max > 0 & bar_min < 0) {
        myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
        myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                      seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
        p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
        
      } else if (bar_max > 0 & bar_min == 0) {
        myColor <- colorRampPalette(c("white", "red"))(paletteLength)
        p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      } else {
        myColor <- colorRampPalette(c("blue","white"))(paletteLength)
        p <- pheatmap(p.yes.rr, color=myColor, fontsize = pheatmap_fontsize, 
                      treeheight_row = treeheight, treeheight_col = treeheight)
      }
      
    } else {
      paletteLength <- 50
      myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
      myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                    seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
      
      p <- pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                    treeheight_row = treeheight, treeheight_col = treeheight)
    }
    
  } else {
    p.yes.rr <- as.matrix(p.yes.rr)
    p = corrplot(p.yes.rr, type="upper", order="hclust", insig = "blank",tl.col = "black",na.label = "o")
  }
  
  detach("package:igraph", unload = TRUE, force = T)
  return(c(list(p = p), list(cor_matrix = p.yes.rr), list(cor_matrix2 = correlation_data_all_2)))
}




