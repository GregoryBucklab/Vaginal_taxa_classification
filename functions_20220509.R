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


######### diffential abundance2 ### no plot, output csv no matter change or not ####################
dif_abundance2 <- function(reads_table,metadata) {
  
  conds <- metadata
  
  x <- aldex.clr(reads_table, conds, mc.samples=128, denom="all", verbose=F)
  
  # paired Wilcoxon Rank Sum test and Welch's t-test
  x.tt <- aldex.ttest(x, paired.test= F)
  
  x.effect <- aldex.effect(x)
  
  x.all <- data.frame(cbind(x.tt,x.effect))
  
  return(list(data = x.all))
  
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




######### RNA-seq differencial analysis #########
Deseq2 <- function(cts, design, treeheight = 50, pheatmap_fontsize = 5) {
  
  
  factor_number <- ncol(design)
  
  if (factor_number == 1) {
    ################ DESeq2 ##################
    # change column names
    col_name <- colnames(design)
    factor_1 = col_name[1]
    
    colnames(design) <- c("factor_1")
    
    design$factor_1 <- as.factor(design$factor_1)
    levels(design$factor_1)  
    
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = design,
                                  design= ~ factor_1)  
    ref1 = toString(design$factor_1[1])
    
    dds$factor_1 <- relevel(dds$factor_1, ref = ref1)   
    
    dds <- DESeq(dds)

    res=resultsNames(dds) # lists the coefficients
    
    res1 <- gsub("factor_1", factor_1, res)
    res1
    
    n=length(res)
    
    lable_list = as.matrix(read.csv('lable.csv', header = F))
    for (i in 2:n)   {
      res2 <- results(dds, name=res[i])
      write.csv(res2,file = paste(res1[i],".csv"))
      
      # volcano plot
      if (i != 4) {
        volvano_plot <- as.data.frame(res2@listData)
        
        volvano_plot$gene_id = row.names(cts)
        volvano_plot$lable_name = NA
        volvano_plot$lable_name[volvano_plot$gene_id %in% lable_list] <- volvano_plot$gene_id[volvano_plot$gene_id %in% lable_list]
        
        volvano_plot$color <- 'grey'
        volvano_plot$color[volvano_plot$log2FoldChange >= 1 & volvano_plot$padj <= 0.05] = 'green'
        volvano_plot$color[volvano_plot$log2FoldChange <= -1 & volvano_plot$padj <= 0.05] = 'red'
        
        
        volvano_plot$padj <- -log10(volvano_plot$padj)
        
        p1=ggplot(data=volvano_plot, aes(x=log2FoldChange, y=padj, color = color, label=lable_name)) + 
          geom_point() + 
          scale_color_manual(values=c("green", "grey", "red")) +
          theme(legend.position = "none") +
          geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
        
        x_limit_neg <- min(quantile(volvano_plot$log2FoldChange,.01, na.rm =T),-2)
        x_limit_pos <- max(quantile(volvano_plot$log2FoldChange,.99, na.rm =T),2)
        volvano_plot$log2FoldChange[volvano_plot$log2FoldChange < x_limit_neg] = x_limit_neg
        volvano_plot$log2FoldChange[volvano_plot$log2FoldChange > x_limit_pos] = x_limit_pos
        
        y_limit_pos <- max(quantile(volvano_plot$padj,.99, na.rm =T),2)
        volvano_plot$padj[volvano_plot$padj > y_limit_pos] = y_limit_pos
        
        p2=ggplot(data=volvano_plot, aes(x=log2FoldChange, y=padj, color = color, label=lable_name)) + 
          geom_point() + 
          scale_color_manual(values=c("green", "grey", "red")) +
          theme(legend.position = "none") +
          geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
        
      }
    }
    
    ############### Data transformations and visualization ###############
    # transfer data using a variance stabilizing transformation
    vsd <- varianceStabilizingTransformation(dds)  
    
    # this gives log2(n + 1)
    ntd <- normTransform(dds)
    
    ############### Heatmap of the count matrix ###############
    select <- order(rowMeans(counts(dds,normalized=TRUE)),
                    decreasing=TRUE)[1:500]
    df <- as.data.frame(colData(dds)[1])
    colnames(df) <- c(factor_1)
    pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE, annotation_col=df, 
             treeheight_row = treeheight, treeheight_col = treeheight)
    
    select <- order(rowMeans(counts(dds,normalized=TRUE)),
                    decreasing=TRUE)[1:nrow(cts)]
    df <- as.data.frame(colData(dds)[1])
    colnames(df) <- c(factor_1)
    pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE, annotation_col=df, fontsize = pheatmap_fontsize, 
             treeheight_row = treeheight, treeheight_col = treeheight)
    
    
    
    ############### Heatmap of the sample-to-sample distances ###############
    sampleDists <- dist(t(assay(vsd)))
    
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$factor_1)     
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors, 
             treeheight_row = treeheight, treeheight_col = treeheight)
    
    
    ############### pincipal component plot of the samples ###############
    pcaData <- plotPCA(vsd, intgroup=c("factor_1"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    #pcaData <- avgdist(cts, sample = )
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    p3=ggplot(pcaData, aes(PC1, PC2, color=factor_1)) +
      geom_point(size=3) +
      labs(color=factor_1) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()
    
    
    
  } else if (factor_number == 2) {
    # change column names
    col_name <- colnames(design)
    factor_1 = col_name[1]
    factor_2 = col_name[2]
    
    colnames(design) <- c("factor_1","factor_2")
    
    design$factor_1 <- as.factor(design$factor_1)
    design$factor_2 <- as.factor(design$factor_2)
    levels(design$factor_1)                                       
    levels(design$factor_2)                                           
    
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = design,
                                  design= ~ factor_1 + factor_2 + factor_2:factor_1 )  
    ref1 = toString(design$factor_1[1])
    ref2 = toString(design$factor_2[1])
    
    dds$factor_1 <- relevel(dds$factor_1, ref = ref1)   
    dds$factor_2 <- relevel(dds$factor_2, ref = ref2)   
    
    dds <- DESeq(dds)
    
    
    res=resultsNames(dds) # lists the coefficients
    
    res1 <- gsub("factor_1", factor_1, res)
    res1 <- gsub("factor_2", factor_2, res1)
    res1
    
    n=length(res)
    
    lable_list = as.matrix(read.csv('lable.csv', header = F))
    for (i in 2:n)   {
      res2 <- results(dds, name=res[i])
      write.csv(res2,file = paste(res1[i],".csv"))
      
      # volcano plot
      if (i != 4) {
        volvano_plot <- as.data.frame(res2@listData)
        
        volvano_plot$gene_id = row.names(cts)
        volvano_plot$lable_name = NA
        volvano_plot$lable_name[volvano_plot$gene_id %in% lable_list] <- volvano_plot$gene_id[volvano_plot$gene_id %in% lable_list]
        
        volvano_plot$color <- 'grey'
        volvano_plot$color[volvano_plot$log2FoldChange >= 1 & volvano_plot$padj <= 0.05] = 'green'
        volvano_plot$color[volvano_plot$log2FoldChange <= -1 & volvano_plot$padj <= 0.05] = 'red'
        
        
        volvano_plot$padj <- -log10(volvano_plot$padj)
        
        p4=ggplot(data=volvano_plot, aes(x=log2FoldChange, y=padj, color = color, label=lable_name)) + 
          geom_point() + 
          scale_color_manual(values=c("green", "grey", "red")) +
          theme(legend.position = "none") +
          geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 
        
        x_limit_neg <- min(quantile(volvano_plot$log2FoldChange,.01, na.rm =T),-2)
        x_limit_pos <- max(quantile(volvano_plot$log2FoldChange,.99, na.rm =T),2)
        volvano_plot$log2FoldChange[volvano_plot$log2FoldChange < x_limit_neg] = x_limit_neg
        volvano_plot$log2FoldChange[volvano_plot$log2FoldChange > x_limit_pos] = x_limit_pos
        
        y_limit_pos <- max(quantile(volvano_plot$padj,.99, na.rm =T),2)
        volvano_plot$padj[volvano_plot$padj > y_limit_pos] = y_limit_pos
        
        p5=ggplot(data=volvano_plot, aes(x=log2FoldChange, y=padj, color = color, label=lable_name)) + 
          geom_point() + 
          scale_color_manual(values=c("green", "grey", "red")) +
          theme(legend.position = "none") +
          geom_text_repel(aes(label = lable_name),size = 3, colour = 'black') 

      }
    }
    
    vsd <- varianceStabilizingTransformation(dds)  # transfer data using a variance stabilizing transformation
    
    ntd <- normTransform(dds)
    
    select <- order(rowMeans(counts(dds,normalized=TRUE)),
                    decreasing=TRUE)[1:20]
    df <- as.data.frame(colData(dds)[,c("factor_1","factor_2")])
    colnames(df) <- c(factor_1,factor_2)
    pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE, annotation_col=df, 
             treeheight_row = treeheight, treeheight_col = treeheight)
    
    select <- order(rowMeans(counts(dds,normalized=TRUE)),
                    decreasing=TRUE)[1:nrow(cts)]
    df <- as.data.frame(colData(dds)[,c("factor_1","factor_2")])
    colnames(df) <- c(factor_1,factor_2)
    pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
             cluster_cols=TRUE, annotation_col=df, 
             treeheight_row = treeheight, treeheight_col = treeheight)
    
    
    sampleDists <- dist(t(assay(vsd)))
    
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$factor_1, vsd$factor_2, sep="_")     
    colnames(sampleDistMatrix) <- NULL
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows=sampleDists,
             clustering_distance_cols=sampleDists,
             col=colors, 
             treeheight_row = treeheight, treeheight_col = treeheight)
    
    
    pcaData <- plotPCA(vsd, intgroup=c("factor_2", "factor_1"), returnData=TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    
    p6=ggplot(pcaData, aes(PC1, PC2, color=factor_2, shape=factor_1)) +
      geom_point(size=3) +
      labs(shape=factor_1, color=factor_2) +
      xlab(paste0("PC1: ",percentVar[1],"% variance")) +
      ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
      coord_fixed()
    
    # output
    output = c(p1 = list(p1),p2 = list(p2),p3 = list(p3),p4 = list(p4),p5 = list(p5),p6 = list(p6))
    return(output)
    
  } else {
    print("Only one or two factors are avaiable in the script")
    return(NA)
  }
}

  
######### diversity analysis #########
diversity_analysis <- function(reads_table_ori,metadata_ori, key_words) {
  # alpha diversity 
  reads_table = t(reads_table_ori)
  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  alpha <- data.frame(diversity(reads_table))
  alpha$Evenness <- alpha$diversity.reads_table./log(specnumber(reads_table))
  alpha$Ovserved_phage <- colSums(t(reads_table) != 0)
  colnames(alpha) = c('Shannon','Evenness','Ovserved_taxa')
  
  alpha$Group = metadata
  
  ggplot(data=alpha, aes(x=Group, y=Shannon)) +
    geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
    geom_jitter(shape=16, size = 0.5)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme_classic()
  ggsave(paste0(key_words,'_Shannon.pdf'),width=2, height=3)
  
  pvalue_output = matrix(data = NA, ncol =2, nrow =0)
  
  pvalue_output_2 = c(paste0(key_words,'_Shannon'), kruskal.test(Shannon~Group, alpha)$p.value)
  pvalue_output = rbind(pvalue_output,pvalue_output_2)
  
  ggplot(data=alpha, aes(x=Group, y=Evenness)) +
    geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
    geom_jitter(shape=16, size = 0.5)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme_classic()
  ggsave(paste0(key_words,'_Evenness.pdf'),width=2, height=3)
  
  pvalue_output_2 = c(paste0(key_words,'_Evenness'), kruskal.test(Evenness~Group, alpha)$p.value)
  pvalue_output = rbind(pvalue_output,pvalue_output_2)
  
  ggplot(data=alpha, aes(x=Group, y=Ovserved_taxa)) +
    geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
    geom_jitter(shape=16, size = 0.5)+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    theme_classic()
  ggsave(paste0(key_words,'_Ovserved_taxa.pdf'),width=2, height=3)
  
  pvalue_output_2 = c(paste0(key_words,'_Ovserved_taxa'), kruskal.test(Ovserved_taxa~Group, alpha)$p.value)
  pvalue_output = rbind(pvalue_output,pvalue_output_2)
  
  
  
  
  # beta diversity
  beta = beta_diversity(reads_table_ori, metadata = metadata_ori, factor_name = 'Group', order = NA, NMDS_skip = F, ref_group = NA, rarefy_to = NA, pheatmap_fontsize = 5,treeheight = 50)
  beta$within_dis_p+
    theme_classic()
  ggsave(paste0(key_words,'_beta_within_dis.pdf'),width=3, height=3)
  
  pvalue_output_2 = c(paste0(key_words,'_within_dis_sig'), beta$within_dis_sig[1,2])
  pvalue_output = rbind(pvalue_output,pvalue_output_2)
  pvalue_output_2 = c(paste0(key_words,'_group_dis_sig'), beta$group_dis_sig[1,2])
  pvalue_output = rbind(pvalue_output,pvalue_output_2)
  
  beta$NMDS
  ggsave(paste0(key_words,'_NMDS.pdf'),width=4, height=3)
  
  beta$NMDS_2
  ggsave(paste0(key_words,'_NMDS_2.pdf'),width=4, height=3)
  
  metadata = as.data.frame(metadata_ori)
  pvalue <- adonis2(reads_table ~ ., data = metadata, method = "bray")  
  pvalue[1,5]
  pvalue_output_2 = c(paste0(key_words,'_Adonis2'), pvalue[1,5])
  pvalue_output = rbind(pvalue_output,pvalue_output_2)
  
  tsne <- Rtsne(reads_table, dims = 2, perplexity=50, verbose=TRUE, max_iter = 5000)
  pic <- tsne$Y
  pic <- data.frame(pic,metadata$metadata_ori)
  
  colnames(pic) <- c('X1','X2','Factor')
  
  ggplot(pic, aes(X1, X2,color = Factor))  +
    geom_point(size=0.5) +
    xlab(paste0("tsne1")) +
    ylab(paste0("tsne2"))+
    stat_ellipse(type = "t") + 
    coord_fixed()+ 
    theme(
      axis.title.x = element_text( size=7),
      axis.title.y = element_text( size=7),
      legend.text = element_text(size=7),
      legend.title = element_text(size=7),
      plot.title = element_text(hjust = 0.5, size = 7)
    ) 
  ggsave(paste0(key_words,'_tsne.pdf'), width=4, height=4)
  
  write.csv(pvalue_output,paste0(key_words,'_pvalues.csv'))
  
  # abundance
  abundance = dif_abundance(reads_table_ori, metadata_ori, 
                             pvalue_th = 0.05, fold_change_th = 1, paired_test = F, style = 1)
  
  if (!is.na(abundance)) {
    abundance$p
    ggsave(paste0(key_words,'_abundance.pdf'), width=4, height=4)

    write.csv(abundance$data,paste0(key_words,'_abundance.csv'))
    
  }
  
  
}



######### ML models ######
#input = metadata; output= icu_metadata
FilterFeatures <- function(input, output)      # remove columns with Standard deviation == 0 
{
  # remove sds = 0
  library(matrixStats)
  
  # remove NA more than 1/3
  keep = colSums(is.na(input)) <= nrow(input)/3
  input = input[,keep]
  
  # remove input with less than 2 levels
  keep = sapply(1:ncol(input), function(j) (length(unique(input[,j])) >= 2 ))
  input = input[,keep]
  
  # remove factor not significantly changes 
  #  keep = sapply(1:ncol(input), function(j) ( kruskal.test(input[,j] ~ as.factor(output))$p.value <=0.05))
  keep = sapply(1:ncol(input), function(j) ( cor.test(input[,j], as.numeric(as.factor(output)), method = 'spearman')$p.value <=0.05))
  
  keep = keep & !is.na(keep)
  input = input[,keep]
  
  return(input)
}

# input num; output binary
xxx<-function(X, Ccohort, i, ntree=500) # X is the input data; Ccohort is the output data; i is seq(length(alpha.shannon_diversity))
{
  library(randomForest)
  library(forestError)
  set.seed(2018+123*i)
  
  if (sum(is.na(X)) >0) {
    X = apply(X, 2, as.numeric)
    Ccohort = as.factor(Ccohort)
    data = cbind(X,Ccohort)
    data_2 <- rfImpute(Ccohort ~ ., data)    # Impute missing values in predictor data using proximity from randomForest.
    X = data_2[,-1]
    X = as.matrix(X)
  } 
  
  XX=X[-i,]    # XX input data after removal of sample i
  CC=Ccohort[-i]    # CC Ccohort data after removal of sample i
  XT=rbind(X[i,], X[i,])   # XT sample i input data  two repeats 
  pvs2=sapply(1:ncol(X), function(j) ( kruskal.test( XX[,j] ~ as.factor(CC))$p.value ) )   # !!! Ccohort is not an continuous dataset, hence it should be converted to factors   # compare input data with cohort 1 vs. 0 in every column 
  cutoff=0.01
  ColFeatInd2=which(pvs2<cutoff)  # find input data columns with significant difference in two cohorts
  while(length(ColFeatInd2)<2)
  {
    ColFeatInd2=which(pvs2<cutoff)
    cutoff=cutoff+0.05
  }
  
  RF2=randomForest(XX[, ColFeatInd2], as.factor(CC), ntree = ntree, keep.inbag = TRUE)  # !!! establish randomForest model. randomForest implements Breiman's random forest algorithm (based on Breiman and Cutler's original Fortran code) for classification and regression. Input is all different input columns and output is two cohorts  
  
  err = quantForestError(RF2, XX[, ColFeatInd2], XT[, ColFeatInd2], ) # conditional misclassification rate
  err = err$mcr[1]
  
  ret=list()
  ret$RF = RF2
  ret$err = err
  ret$p2 = predict(RF2, XT[, ColFeatInd2], type='prob')[1]  # !!! predict using RF2 randomForest model, input is sample i. output is expected value.
  ret$ptrain2 = predict(RF2, XX[, ColFeatInd2], type='prob')   # predict using RF2 randomForest model, input are samples except for i.
  ret$coef2=rep(0, ncol(XX))
  ret$coef2[ColFeatInd2]=importance(RF2)   # importance of used input columns in RF2 model
  return(ret)
}

twoLayer<-function(X, Ccohort, i, ntree=500)
{
  ret=list()
  ret$Layer1=xxx(X, Ccohort, i)   # get RF2 randomForest model and prediction of sample i and other samples
  prdV=list()
  for(j in seq(nrow(X[-i,])))   # 1:(all sample number - 1) loop 
    prdV[[j]]=xxx(X[-i,], Ccohort[-i], j)  #  
  ret$Layer2=prdV
  return(ret)
}

# xxx and get output 
xxx_get_result <- function(input, output, ntree = 500) {
  
  xxx<-function(X, Ccohort, i, ntree=500) # X is the input data; Ccohort is the output data; i is seq(length(alpha.shannon_diversity))
  {
    library(randomForest)
    library(forestError)
    set.seed(2018+123*i)
    
    if (sum(is.na(X)) >0) {
      X = apply(X, 2, as.numeric)
      Ccohort = as.factor(Ccohort)
      data = cbind(X,Ccohort)
      data_2 <- rfImpute(Ccohort ~ ., data)    # Impute missing values in predictor data using proximity from randomForest.
      X = data_2[,-1]
      X = as.matrix(X)
    } 
    
    XX=X[-i,]    # XX input data after removal of sample i
    CC=Ccohort[-i]    # CC Ccohort data after removal of sample i
    XT=rbind(X[i,], X[i,])   # XT sample i input data  two repeats 
    pvs2=sapply(1:ncol(X), function(j) ( kruskal.test( XX[,j] ~ as.factor(CC))$p.value ) )   # !!! Ccohort is not an continuous dataset, hence it should be converted to factors   # compare input data with cohort 1 vs. 0 in every column 
    cutoff=0.01
    ColFeatInd2=which(pvs2<cutoff)  # find input data columns with significant difference in two cohorts
    while(length(ColFeatInd2)<2)
    {
      ColFeatInd2=which(pvs2<cutoff)
      cutoff=cutoff+0.05
    }
    
    RF2=randomForest(XX[, ColFeatInd2], as.factor(CC), ntree = ntree, keep.inbag = TRUE)  # !!! establish randomForest model. randomForest implements Breiman's random forest algorithm (based on Breiman and Cutler's original Fortran code) for classification and regression. Input is all different input columns and output is two cohorts  
    
    err = quantForestError(RF2, XX[, ColFeatInd2], XT[, ColFeatInd2], ) # conditional misclassification rate
    err = err$mcr[1]
    
    ret=list()
    ret$RF = RF2
    ret$err = err
    ret$p2 = predict(RF2, XT[, ColFeatInd2], type='prob')[1]  # !!! predict using RF2 randomForest model, input is sample i. output is expected value.
    ret$ptrain2 = predict(RF2, XX[, ColFeatInd2], type='prob')   # predict using RF2 randomForest model, input are samples except for i.
    ret$coef2=rep(0, ncol(XX))
    ret$coef2[ColFeatInd2]=importance(RF2)   # importance of used input columns in RF2 model
    return(ret)
  }
  
  prd=foreach(i=seq(nrow(input))) %dopar% xxx(input, output, i, ntree= ntree)   # foreach(%dopar%) %dopar% evaluates it in parallel, while %do% evaluates the expression sequentially. 
  
  ppp=vector()
  ML_error = vector()
  for(i in seq(nrow(input))) {
    ppp[i] = prd[[i]]$p2
    ML_error[i] = prd[[i]]$err
  }
  err_median = median(ML_error)
  
  roc = roc(output, ppp)
  data = cbind(roc$sensitivities,roc$specificities)
  data = as.data.frame(data)
  colnames(data) = c('Sensitivities','Specificities')
  data = data[order(data$Specificities, decreasing = T),]
  p_roc = ggplot(data, aes(x = Sensitivities, y = Specificities)) +
    geom_point(size = 0.3)+geom_step()+
    theme(axis.title = element_text(size = 6), 
          axis.text = element_text(size = 6), 
          legend.text = element_text(size = 6), 
          legend.title = element_text(size = 6)) + theme_bw() 
  
  auc = roc(output, ppp)$auc      # It builds a ROC curve and returns a “roc” object, a list of class “roc”.
  pvalue = -log10(wilcox.test(ppp ~ as.factor(output))$p.value)   # test correlation between expected values and true values (PTB yes or no).
  
  importance_ML = as.data.frame(matrix(data = NA, ncol = length(prd), nrow = ncol(input)))
  for (a in seq(length(prd))) {
    importance_ML[,a] = prd[[a]]$coef2
  }
  importance_ML = rowMeans(importance_ML)
  
  importance_ML <- data.frame(name=colnames(input),
                              MeanDecreaseGini = importance_ML)
  
  importance_ML$Association = NA
  importance_ML$Yes_quantile = NA
  importance_ML$No_quantile = NA
  
  for (b in seq(nrow(importance_ML))) {
    data_2 = input[,b]
    data_3 = data_2[output == 'Yes']
    data_4 = data_2[output == 'No']
    
    quantile1 = quantile(data_3[!is.na(data_3)])
    quantile2 = quantile(data_4[!is.na(data_4)])
    
    if (quantile1[3] > quantile2[3]) {
      importance_ML$Association[b]= 'Positive'
    } else {
      importance_ML$Association[b]= 'Negative'
    }
    importance_ML$Yes_quantile[b] = paste(as.character(quantile1), collapse = '; ')
    importance_ML$No_quantile[b] = paste(as.character(quantile2), collapse = '; ')
  }
  importance_ML = importance_ML[order(importance_ML$MeanDecreaseGini, decreasing = T),]
  
  ML_result =list()
  ML_result$prd = prd
  ML_result$p_roc =p_roc
  ML_result$auc=auc
  ML_result$err_median = err_median
  ML_result$pvalue=pvalue
  ML_result$importance_ML=importance_ML
  return(ML_result)
  
}

# input num/cha; output interval
xxx_interval <-function(X, Ccohort, type = 'num', mtry=5, min.node.size = 5) # X is the input data; Ccohort is the output data; i is seq(length(alpha.shannon_diversity))
{
  
  library(caret)
  
  set.seed(2022)
  
  if (sum(is.na(X)) >0) {
    if (type == 'num') {
      X = apply(X, 2, as.numeric)
      data = cbind(X,Ccohort)
      
    } else {
      X = lapply( X, factor)
      X = as.data.frame(X)
      data = sapply(1:ncol(X), function(j) ( as.numeric(X[,j]) ))
      data = as.data.frame(data)
      colnames(data) = colnames(X)
    }
    
    
    if (type == 'num') {
      library(randomForest)
      data_2 <- rfImpute(Ccohort ~ ., data)    # Impute missing values in predictor data using proximity from randomForest.
      X = data_2[,-1]
    } else {
      # With mice
      library(mice)
      #      library(missMDA)
      x.impmi<- mice(data, m = 2, printFlag = FALSE)
      x.impmi_2 = x.impmi$imp
      
      x.impmi_name = names(data)  
      for (b in 1:length(x.impmi_2)) {
        x.impmi_3 = x.impmi_2[[b]]
        colnames(x.impmi_3) = c('A','B')
        x.impmi_3$result = (x.impmi_3$A + x.impmi_3$B)/2
        
        n = which(colnames(data) == x.impmi_name[b])
        
        for (c in 1: nrow(x.impmi_3)) {
          m = as.numeric(row.names(x.impmi_3)[c])
          data[m,n]= x.impmi_3$result[c]
        }
        
      }
      X = data
    }
  } 
  
  # we are not going to do any cross-validatin and rely on OOB error
  trctrl <- trainControl(method = "none")
  data = cbind(Ccohort,X)
  # we will now train random forest model
  rfregFit <- train(Ccohort~., data = data, method = "ranger",importance="permutation",
                    tuneGrid = data.frame(mtry=mtry,
                                          min.node.size = min.node.size,
                                          splitrule="variance"))
  return(rfregFit)
}

##### glm ###
#X=input; Ccohort = output; i=2
xxx_glm <-function(X, Ccohort, i) # X is the input data; Ccohort is the output data; i is seq(length(alpha.shannon_diversity))
{
  
  XX=X[-i,]    # XX input data after removal of sample i
  CC=Ccohort[-i]    # CC Ccohort data after removal of sample i
  XT=rbind(X[i,], X[i,])   # XT sample i input data  two repeats 
  
  data = cbind(XX, CC)
  RF2= glm(CC~., data = data, method = "glm.fit", family = gaussian)
  
  ret=list()
  ret$p2 = predict.glm(RF2, XT)[1]  # !!! predict using RF2 randomForest model, input is sample i. output is expected value.
  ret$ptrain2 = predict(RF2, XX)   # predict using RF2 randomForest model, input are samples except for i.
  ret$coef2= RF2$coefficients
  return(ret)
}

##### gnm ###
xxx_gnm <-function(X, Ccohort, i) # X is the input data; Ccohort is the output data; i is seq(length(alpha.shannon_diversity))
{
  library(gnm)
  XX=X[-i,]    # XX input data after removal of sample i
  CC=Ccohort[-i]    # CC Ccohort data after removal of sample i
  XT=rbind(X[i,], X[i,])   # XT sample i input data  two repeats 
  
  data = cbind(XX, CC)
  RF2= gnm(CC~., data = data, method = "gnmFit", family = gaussian)
  
  ret=list()
  ret$p2 = predict(RF2, XT)[1]  # !!! predict using RF2 randomForest model, input is sample i. output is expected value.
  ret$ptrain2 = predict(RF2, XX)   # predict using RF2 randomForest model, input are samples except for i.
  ret$coef2= RF2$coefficients
  return(ret)
}


######### sample size power #####

# fit distribution
#library(fitdistrplus)
# descdist(data_1, discrete=F)

#library(pwr)
# know power and calculate sample size
# d = (mean(data_1)-mean(data_2))/sqrt(((sd(data_1))^2 + (sd(data_2))^2)/2) # Effect size
# pwr.t.test(d=d, sig.level=0.05, power=0.80, type= "two.sample", alternative="greater")

#library(wmwpow)
# know sample size and calculate power
# shiehpow(n = length(data_1), m = length(data_1), p = 0.80, alpha = 0.05, dist = "norm", sides = "two.sided")
