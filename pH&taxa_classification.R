##### get and prepare data, 'BT' are not involved #########
setwd('/Users/binzhu/secure/vamp/home/metadata/vahmp/stirrups')
data <- read.table('2017-10-10-VT_Uclust_filtered.txt', sep = '\t')
setwd('/Users/binzhu/secure/vamp/home/metadata/vahmp/healthHistory')
metadata <- read.table('2017-03-24vampClinicalData.txt', sep = '\t') 
setwd('/Users/binzhu/Desktop/pH/results')
data <- data[-1,]

### reads table ###
sample_name <- unique(data$V1)
species_name <- unique(data$V2)
data$V2 <- as.character(data$V2)
data$V1 <- as.character(data$V1)
data$V4 <- as.numeric(as.character(data$V4))
keep <- data$V3 == 'AT'
data <- data[keep,]

reads_table <- matrix(0, ncol = length(sample_name), nrow = length(species_name))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = species_name 
colnames(reads_table) = sample_name

for (a in 1: nrow(data)) {
   
   column_num <- which(sample_name == data[a,1])
   row_num <- which(species_name == data[a,2])
   
   reads_table[row_num,column_num] =  data$V4[a]
   
}
reads_table <- as.data.frame(reads_table)
row.names(reads_table)[row.names(reads_table) == "TM7_OTU-H1"] = "BVAB_TM7"

# change to abundance
reads_table_all = reads_table

reads_table_all_abundance <- as.data.frame(matrix(data = 0, nrow = nrow(reads_table_all), ncol =ncol(reads_table_all)))
row.names(reads_table_all_abundance) <- row.names(reads_table_all)
colnames(reads_table_all_abundance) <- colnames(reads_table_all)

max_species <- as.data.frame(matrix(data = NA, nrow = ncol(reads_table_all), ncol =1))
sum_number <- as.data.frame(matrix(data = NA, nrow = ncol(reads_table_all), ncol =1))

for (a in 1:dim(reads_table_all)[2]) {
   sum_number[a,1] = sum(reads_table_all[,a])
   reads_table_all_abundance[,a] <- reads_table_all[,a] / sum_number[a,1]
   max_species[a,1] <- row.names(reads_table_all_abundance)[which.max(reads_table_all_abundance[,a])]
   
}

Top_abundant_taxa = rowSums(reads_table_all_abundance)/ncol(reads_table_all_abundance)
Top_abundant_taxa = as.data.frame(sort(Top_abundant_taxa, decreasing = T))
Top_abundant_taxa$taxa = row.names(Top_abundant_taxa)

# Species with relative abundance >= 1% and species with the same genus as top abundant species kept
# Other species combine to genus level
Top_abundant_taxa = Top_abundant_taxa[Top_abundant_taxa$`sort(Top_abundant_taxa, decreasing = T)` >= 0.01,]
Top_abundant_taxa_genus = str_replace_all(Top_abundant_taxa$taxa, '_.*', '')
Top_abundant_taxa_genus = unique(Top_abundant_taxa_genus)

keep = str_detect(row.names(reads_table_all), 'Lactobacillus|Gardnerella|Lachnospiraceae|Atopobium|Sneathia|Prevotella')
row.names(reads_table_all)[keep]

reads_table_1 = reads_table_all[keep,]
reads_table_2 = reads_table_all[!keep,]

keep = row.names(reads_table_1) %in% Top_abundant_taxa$taxa
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

# species threshold
reads_table = prepare_reads_table_2(reads_table, total_reads_threshold = 5000, species_threshold = 0.001, mc.cores = 8) 
reads_table = reads_table$reads_table
row.names(reads_table)[row.names(reads_table) == 'Lactobacillus_crispatus'] = 'Lactobacillus_crispatus'
row.names(reads_table)[row.names(reads_table) == 'Lactobacillus_gasseri'] = 'Lactobacillus_gasseri'
row.names(reads_table)[row.names(reads_table) == 'Ca._L._vaginae'] = 'Ca._L._vaginae'
row.names(reads_table)[row.names(reads_table) == 'TM7'] = 'BVAB_TM7'
#row.names(reads_table)[row.names(reads_table) == 'Coriobacteriaceae_spp.'] = 'Coriobacteriaceae'
reads_table_all = reads_table

### prepare metadata ###
# metadata
metadata$V1 = str_replace_all(metadata$V1, 'S', 'VT')

sample_name <- unique(metadata$V1)
survey <- unique(metadata$V2)

reads_table <- matrix(NA, nrow = length(sample_name), ncol = length(survey))     # create reads table   sample name as columns and species name as rows
row.names(reads_table) = sample_name 
colnames(reads_table) =survey 

metadata$V4 <- as.character(metadata$V4)
metadata$V1 <- as.character(metadata$V1)
metadata$V2 <- as.character(metadata$V2)

for (a in 1: nrow(metadata)) {
   kitID_num <- which(sample_name == metadata[a,1])
   survey_num <- which(survey == metadata[a,2])
   reads_table[kitID_num,survey_num] = metadata[a,4]
}

metadata_all = reads_table
metadata_all = as.data.frame(metadata_all)

# match reads table with metadata
keep = colnames(reads_table_all) %in% row.names(metadata_all)
sum(keep)
reads_table_all = reads_table_all[,keep]
keep = row.names(metadata_all) %in% colnames(reads_table_all)
sum(keep)
metadata_all = metadata_all[keep,]
reads_table_all = reads_table_all[row.names(metadata_all)]

rm(reads_table, metadata,cts)

### modify metadata
metadata_all = metadata_all[,!str_detect(colnames(metadata_all),'reason|child|dob')]
metadata_all[metadata_all == 'NS'] = NA
metadata_all[metadata_all == 'N/A'] = NA
# alcohol_30days
metadata_all[metadata_all == '0'] = 0
metadata_all[metadata_all == '1-2'] = 1
metadata_all[metadata_all == '3'] = 2
metadata_all[metadata_all == '6-9'] = 3
metadata_all[metadata_all == '10-19'] = 4
metadata_all[metadata_all == '20+'] = 5
metadata_all[metadata_all == 'No Answer'] = NA
# number_sexpartners_life
metadata_all[metadata_all == '0'] = 0
metadata_all[metadata_all == '1'] = 1
metadata_all[metadata_all == '2'] = 2
metadata_all[metadata_all == '3-5'] = 3
metadata_all[metadata_all == '6-10'] = 4
metadata_all[metadata_all == '11-20'] = 5
metadata_all[metadata_all == '21+'] = 6
# education
metadata_all[metadata_all == 'two year college'] = 3
metadata_all[metadata_all == 'doctoral degree'] = 6
metadata_all[metadata_all == 'high school'] = 2
metadata_all[metadata_all == 'less than high school'] = 1
metadata_all[metadata_all == 'four year college'] = 4
metadata_all[metadata_all == 'some college'] = 3
metadata_all[metadata_all == 'masters degree'] = 5
# income
metadata_all[metadata_all == '15k to 20k'] = 2
metadata_all[metadata_all == '20k to 40k'] = 3
metadata_all[metadata_all == '40k to 60k'] = 4
metadata_all[metadata_all == '60k to 80k'] = 5
metadata_all[metadata_all == 'more than 80k'] = 6
metadata_all[metadata_all == 'less than 15k'] = 1
# douche_frequency
metadata_all[metadata_all == 'one to three per month'] = 2
metadata_all[metadata_all == 'two to six per week'] = 4
metadata_all[metadata_all == 'less than once per month'] = 1
metadata_all[metadata_all == 'Never'] = 0
metadata_all[metadata_all == 'never'] = 0
metadata_all[metadata_all == 'once per day'] = 5
metadata_all[metadata_all == 'once per week'] = 3
# douche
metadata_all[metadata_all == 'one to three weeks'] = 2
metadata_all[metadata_all == 'one to two days'] = 4
metadata_all[metadata_all == 'more than one month'] = 1
metadata_all[metadata_all == 'more than one year ago or never'] = 0
metadata_all[metadata_all == 'last 24 hours'] = 5
metadata_all[metadata_all == 'three to seven days'] = 3
# bacterial_vaginosis_frequency
metadata_all[metadata_all == 'one'] = 1
metadata_all[metadata_all == 'two to four'] = 2
metadata_all[metadata_all == 'five or more'] = 3
# vigorous_physical
metadata_all[metadata_all == 'zero times'] = 0
metadata_all[metadata_all == 'one to two times'] = 1
metadata_all[metadata_all == 'three to four times'] = 2
metadata_all[metadata_all == 'five to six times'] = 3
metadata_all[metadata_all == 'seven or more times'] = 4
# yogurt
metadata_all[metadata_all == 'less than one per week'] = 1
metadata_all[metadata_all == 'one to two per week'] = 2
metadata_all[metadata_all == 'three to four per week'] = 3
metadata_all[metadata_all == 'one to two per day'] = 4
metadata_all[metadata_all == 'three to four per day'] = 5
metadata_all[metadata_all == 'five plus per day'] = 6
# flow
metadata_all[metadata_all == 'Light'] = 1
metadata_all[metadata_all == 'Medium'] = 2
metadata_all[metadata_all == 'Heavy'] = 3
metadata_all[metadata_all == 'Heavy,Medium,Light'] = NA
metadata_all[metadata_all == 'Heavy,Medium'] = NA
metadata_all[metadata_all == 'Heavy,Light'] = NA
metadata_all[metadata_all == 'Medium,Light'] = NA
# bleedingfreq
metadata_all[metadata_all == 'Less frequently'] = 1
metadata_all[metadata_all == 'Every other cycle'] = 2
metadata_all[metadata_all == 'Every cycle'] = 3
# pillstime
metadata_all[metadata_all == 'More than 2 years'] = 1
metadata_all[metadata_all == '1-2 years'] = 2
metadata_all[metadata_all == '6-12 months'] = 3
metadata_all[metadata_all == '3-6 months'] = 4
metadata_all[metadata_all == '1-3 months'] = 5
# alcohol_30days
metadata_all[metadata_all == 'less than 15 minutes'] = 1
metadata_all[metadata_all == '15 to 30 minutes'] = 2
metadata_all[metadata_all == '30 to 60 minutes'] = 3
metadata_all[metadata_all == 'more than 1 hour'] = 4
# stress
metadata_all[metadata_all == 'much less stress than most people'] = 1
metadata_all[metadata_all == 'a bit less stress than most people'] = 2
metadata_all[metadata_all == 'about average stress'] = 3
metadata_all[metadata_all == 'a bit more stress than most people'] = 4
metadata_all[metadata_all == 'much more stress than most people'] = 5
# twin_meeting
metadata_all[metadata_all == 'less than one per year'] = 1
metadata_all[metadata_all == 'one to five times per year'] = 2
metadata_all[metadata_all == 'six to eleven times per year'] = 3
metadata_all[metadata_all == 'one to three times per month'] = 4
metadata_all[metadata_all == 'one to three times per week'] = 5
metadata_all[metadata_all == 'four or more times per week'] = 6
# worries
metadata_all[metadata_all == 'Much less'] = 1
metadata_all[metadata_all == 'A bit less'] = 2
metadata_all[metadata_all == 'Average'] = 3
metadata_all[metadata_all == 'A bit more'] = 4
metadata_all[metadata_all == 'Much more'] = 5
# something_unexpected_worries
metadata_all[metadata_all == 'Never'] = 0
metadata_all[metadata_all == 'Almost never'] = 1
metadata_all[metadata_all == 'Sometimes'] = 2
metadata_all[metadata_all == 'Fairly often'] = 3
metadata_all[metadata_all == 'Very often'] = 4

reads_table_all = reads_table_all[row.names(metadata_all)]
colnames(reads_table_all) = paste0('Pseudo_sample_id_',seq(ncol(reads_table_all)))     ### run when need local output
row.names(metadata_all) = paste0('Pseudo_sample_id_',seq(ncol(reads_table_all)))     ### run when need local output

c2 = c('0','1','2','3','4','5','6','7','8','9','.','-')
for (a2 in 1:ncol(metadata_all)) {
   b2 = metadata_all[,a2]
   b2 = b2[!is.na(b2)]
   b2 = strsplit(b2,'*')
   b2 = unlist(b2)
   b2 = unique(b2)
   keep = b2 %in% c2
   
   if (sum(!keep) == 0) {
      metadata_all[,a2] = as.numeric(metadata_all[,a2])
   } 
}
x = list()
for (a in seq(ncol(metadata_all))) {
   x[[a]] = unique(metadata_all[,a])
}
rm(reads_table, reads_table_1, reads_table_2,reads_table_3,reads_table_all_abundance, max_species,data, sum_number,x,Top_abundant_taxa)

# prepare data
{
  keep = which(metadata_all$sample_ph == 70)
  metadata_all = metadata_all[-keep,]
  reads_table_all = reads_table_all[,-keep]
  keep = colSums(is.na(metadata_all)) < nrow(metadata_all)/3
  metadata_all = metadata_all[,keep]
}
reads_table_all = reads_table_all[-which(row.names(reads_table_all) == 'Other_Prevotellaceae'),]

write.csv(row.names(reads_table_all),'taxa_list.csv')
write.csv(reads_table_all,'reads_table_all.csv')
write.csv(metadata_all,'metadata_all.csv')









################# diversity #################
reads_table_all_abundance = get_abundance_table(reads_table_all)
type_th = 0.3

# species list and color list
{
   species_list = c("Lactobacillus_crispatus","Lactobacillus_gasseri", 
                    "Lactobacillus_iners","Lactobacillus_jensenii","Atopobium_vaginae","Gardnerella_vaginalis","Ca._L._vaginae", 
                    "Streptococcus", "No_Type","Others","Mycoplasma","Sneathia_amnii", "Prevotella_bivia")

  species_list_2 = c("Lactobacillus crispatus","Lactobacillus gasseri", 
                   "Lactobacillus iners","Lactobacillus jensenii","Atopobium vaginae","Gardnerella vaginalis","Ca. L. vaginae", 
                   "Streptococcus", "No Type","Others","Mycoplasma","Sneathia amnii", "Prevotella bivia")
  
   color_list = c("#FFDC00","#c5b0d5", 
                  "#aec7e8","#49C800","#c49c94","#d62728", "#ff7f0e", 
                  "#7f7f7f","#4d4d4d","#9C9146","#e377c2", "#FF00D1", "#17becf")
   
   color_list_match = c("Lactobacillus crispatus" = "#FFDC00","Lactobacillus jensenii" = "#49C800","Lactobacillus gasseri" = "#c5b0d5", 
                  "Lactobacillus iners"="#aec7e8","Ca. L. vaginae"="#ff7f0e", "Gardnerella vaginalis"="#d62728",  
                  "Sneathia amnii"="#FF00D1", "Prevotella bivia"="#17becf","Atopobium vaginae"="#c49c94",
                  "Mycoplasma"="#e377c2", "Streptococcus"="#7f7f7f","No Type"="#4d4d4d","Others"="#9C9146")
   
}

# get vagitype
{
   reads_table = as.data.frame(t(reads_table_all_abundance))
   reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
   
   mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
   maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
   mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
   mytypes[maxprop < type_th] <- "No_Type"
   uniqueTypes <- unique(mytypes)
   
   data = cbind(colnames(reads_table_all), mytypes)
   write.csv(data,'Vagitype.csv')
   top_vagitype = as.data.frame(table(mytypes))
   top_vagitype = top_vagitype[order(top_vagitype$Freq, decreasing = T),]
   
   write.csv(top_vagitype,'top_vagitype.csv')
   keep = mytypes %in% top_vagitype$mytypes[1:12]
   
   mytypes[!keep] = 'Others'
   metadata_all$Vagitype = mytypes
}

# the vaginal ph associated with vagitypes
{
    reads_table = as.data.frame(t(reads_table_all_abundance))
    reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
    
    mytypes_2 <- apply(reads_table,1,which.max)   # find the most abundant taxa
    maxprop_2 <- reads_table[matrix(c(1:nrow(reads_table),mytypes_2), ncol=2)]  # the abundance of the most abundant taxa
    mytypes_2 <- colnames(reads_table)[mytypes_2]   # find the name of the most abundant taxa
    mytypes_2[maxprop_2 <= 0.3] <- "No_Type"
  
  data = as.data.frame(matrix(data = NA, ncol = 2, nrow = ncol(reads_table_all)))
  colnames(data) = c('Vagitype','Vaginal_pH')
  data$Vagitype = mytypes_2
  data$Vaginal_pH = metadata_all$sample_ph
  unique(data$Vagitype)
  
  data = data[!is.na(data$Vaginal_pH),]
  
  taxa_list = unique(data$Vagitype)
  data_2 = as.data.frame(matrix(data = NA, ncol = 5, nrow = length(taxa_list)))
  colnames(data_2) = c('Predominated_taxa','Median_value','Median_25th_75th_empirical_quartile','Mean_SD','Case_number')
  data_2$Predominated_taxa = taxa_list
  
  for (a in 1:nrow(data_2)) {
    n = which(data$Vagitype == data_2$Predominated_taxa[a])
    data_3 = data$Vaginal_pH[n]
    data_2$Case_number[a] = length(data_3)
    data_3 = quantile(data_3)
    data_2$Median_value[a] = data_3[3]
    data_2$Median_25th_75th_empirical_quartile[a] = paste0(data_3[3], ' (',data_3[2],", ",data_3[4],')')
    data_2$Mean_SD[a] = paste0(round(mean(data$Vaginal_pH[n]), 3), ' (',round(sd(data$Vaginal_pH[n]), 3),')')
    
  }
  write.csv(data_2,'ph_0.3_dominance.csv')
}


# t-SNE
{
   cts2_tsne <- as.data.frame(t(reads_table_all))
   cts2_tsne <- as.data.frame(clr(cts2_tsne))
   tsne <- Rtsne(cts2_tsne, dims = 2, perplexity=500, verbose=TRUE, max_iter = 1500,check_duplicates = FALSE)
   
   pic <- tsne$Y
   pic <- data.frame(pic,mytypes)
   colnames(pic) <- c('X1','X2','Vagitype')
   pic$SampleID <- row.names(cts2_tsne)
   pic$Vagitype <- factor(pic$Vagitype, levels = species_list)
   pic$Vagitype = str_replace_all(pic$Vagitype, "_"," ")
   
   ggplot(pic, aes(X1, X2,color = Vagitype))  +
    geom_point(size=0.5,aes(color = Vagitype)) +
    scale_color_manual(values= color_list_match) +
    coord_flip()+
    xlab("t-SNE1") +
    ylab("t-SNE2") + theme_bw()+
    theme(
         axis.title.x = element_text( size=12),
         axis.title.y = element_text( size=12),
         legend.text = element_text(size=12),
         legend.title = element_text(size=12),
         plot.title = element_text(hjust = 0.5, size = 14)
    )
   ggsave("t_SNE_1.pdf", width=6, height=6)
   
   ggplot(pic, aes(X1, X2,color = Vagitype))  +
     geom_point(size=0.5,aes(color = Vagitype)) +
     scale_color_manual(values= color_list_match) +
     stat_ellipse(type = "euclid",size = 1)+
     coord_flip()+
     xlab("t-SNE1") +
     ylab("t-SNE2") + theme_bw()+
     theme(
       axis.title.x = element_text( size=12),
       axis.title.y = element_text( size=12),
       legend.text = element_text(size=12),
       legend.title = element_text(size=12),
       plot.title = element_text(hjust = 0.5, size = 14)
     )
   ggsave("t_SNE_2.pdf", width=6, height=6)
}


# alpha diversity
{
  reads_table = reads_table_all; metadata = mytypes; factor_name = NA;
  order = NA; NMDS_skip = T; ref_group = NA; rarefy_to = NA; pheatmap_fontsize = 50;treeheight = 10;
  pheatmap_y = T
  
  # rarefy to normalize data
  reads_table = as.data.frame(t(reads_table))    
  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  metadata=as.matrix(metadata)
  
  # alpha
  # calculate diversity
  detach("package:igraph", unload = TRUE, force = T)
  alpha.shannon_diversity <- data.frame(diversity(reads_table))
  alpha.evenness <- alpha.shannon_diversity/log(specnumber(reads_table))
  alpha.ovserved_OTU <- data.frame(colSums(t(reads_table) != 0))
  
  alpha = as.data.frame(matrix(data = NA,ncol=3,nrow = nrow(reads_table)))
  colnames(alpha) = c('alpha.shannon','alpha.evenness','alpha.ovserved_OTU')
  
  alpha$alpha.shannon <- alpha.shannon_diversity$diversity.reads_table.
  alpha$alpha.evenness <- alpha.evenness$diversity.reads_table.
  alpha$alpha.ovserved_OTU <- alpha.ovserved_OTU$colSums.t.reads_table.....0.
  
  metadata_2 = cbind(metadata, alpha)
  
  colnames(metadata_2)[1] = 'Vagitype'
  metadata_2$Vagitype = str_replace_all(metadata_2$Vagitype, "_"," ")
  metadata_2$Vagitype = factor(metadata_2$Vagitype, levels = species_list_2)
  
  ggplot(metadata_2, aes(x=Vagitype, y=alpha.shannon, group = Vagitype, fill = Vagitype))+
    scale_fill_manual(values=color_list)+
    geom_boxplot(outlier.shape=NA)+
    labs(x = NULL, y = "Shannon index", fill=species_list_2)+theme_bw()+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7),
          axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))
  ggsave("shannon_distance.pdf", width=4, height=3)
  
  metadata_2$Vagitype = as.character(metadata_2$Vagitype)
  factor_levels = unique(metadata_2$Vagitype)
  n = length(factor_levels)
  
  Shannon_sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
  colnames(Shannon_sig) = factor_levels
  row.names(Shannon_sig) = factor_levels
  
  for (a in 1:(n-1)) {
    for (b in (a+1) : n) {
      factor_level1 <- subset(metadata_2,  Vagitype == factor_levels[a],
                              drop = TRUE)
      factor_level2 <- subset(metadata_2,  Vagitype == factor_levels[b],
                              drop = TRUE)
      
      Shannon_sig[a,b] <- wilcox.test(factor_level1$alpha.shannon, 
                                      factor_level2$alpha.shannon, paired = F)$p.value

      
    }
  }
  write.csv(Shannon_sig,'Shannon_sig.csv', row.names = T, quote = F)

}

# beta diversity
{
  # prepare data
  {
    reads_table = reads_table_all; metadata = mytypes; factor_name = NA;
    order = NA; NMDS_skip = T; ref_group = NA; rarefy_to = NA; pheatmap_fontsize = 50;treeheight = 10;
    pheatmap_y = T
    
    # rarefy to normalize data
    reads_table = as.data.frame(t(reads_table))    
    reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
    reads_table <- reads_table$otu.tab.rff
    reads_table <- as.data.frame(reads_table)
    
    metadata=as.matrix(metadata)
    
    # Bray_Curtis
    Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
    Bray_Curtis <- as.data.frame(Bray_Curtis)
    
    Bray_Curtis_2 = Bray_Curtis
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
  }
  
  # within group
  {
    keep = group_dis$Source == group_dis$Target
    within_dis = group_dis[keep,]
    keep = within_dis$key != within_dis$key2
    within_dis = within_dis[keep,]
    
    keep = within_dis$Source %in% species_list
    within_dis = within_dis[keep,]
    
    within_dis$Source = factor(within_dis$Source , levels = species_list)
    
    within_dis$Source = str_replace_all(within_dis$Source, "_"," ")
    within_dis$Source = factor(within_dis$Source, levels = species_list_2)
    
    ggplot(within_dis, aes(x=Source, y=value, group = Source, fill = Source))+
      scale_fill_manual(values=color_list)+
      geom_boxplot(outlier.shape=NA)+
      labs(x = NULL, y = "Within group distance", fill=factor_name)+ theme_bw()+
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7),
            axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))
    ggsave("within_group_distance.pdf", width=4, height=3)
    ggsave("within_group_distance_bar_legend.pdf", width=4, height=5)
    
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
    rm(Source, within_dis)
  }
  
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
    
    diag(group_dis_sig) = 1
    pdf("between_group_distance.pdf")
    print(corrplot(distance_median,p.mat = group_dis_sig, method = 'shade', diag = F, type="upper",
                   sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.9,insig = 'label_sig', pch.col = 'grey20', 
                   is.corr=FALSE, col=colorRampPalette(c("red","white"))(100)[c(10:100)], hclust.method = "complete" ,tl.col = "black"))
    dev.off()
    
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
    
    for (a in 1:n) {
      for (b in 1:n) {
        if (a == b) {
          distance_median[a,b] = NA
        }
      }
    }
    
    colnames(distance_median) = str_replace_all(colnames(distance_median),"_"," ")
    row.names(distance_median) = str_replace_all(row.names(distance_median),"_"," ")
    group_dis_2_p = pheatmap(distance_median, cluster_rows=TRUE, show_rownames=TRUE, 
                             cluster_cols=T, show_colnames=T, 
                             color=colorRampPalette(c("red","white"))(100),
                             fontsize = 15, display_numbers = group_dis_sig_2, 
                             treeheight_row = 50, treeheight_col = 50, method = "ward.D",
                             clustering_distance_rows = "euclidean",
                             clustering_distance_cols = "euclidean")
    save_heatmap_pdf(group_dis_2_p, "between_group_distance_2_ward.D_euclidean.pdf", width=7, height=7)
    write.csv(distance_median,'between_group_distance_median.csv', row.names = T, quote = F)
    
    distance_median_2 = distance_median[-12,]
    distance_median_2 = distance_median_2[,-12]
    group_dis_sig_3 = group_dis_sig_2
    group_dis_sig_3 = group_dis_sig_3[-12,]
    group_dis_sig_3 = group_dis_sig_3[,-12]
    group_dis_2_p = pheatmap(distance_median_2, cluster_rows=TRUE, show_rownames=TRUE, 
                             cluster_cols=T, show_colnames=T, 
                             color=colorRampPalette(c("red","white"))(100),
                             fontsize = 15, display_numbers = group_dis_sig_3, 
                             treeheight_row = 50, treeheight_col = 50, method = "ward.D",
                             clustering_distance_rows = "euclidean",
                             clustering_distance_cols = "euclidean")
    save_heatmap_pdf(group_dis_2_p, "between_group_distance_2_ward.D_euclidean_no_others.pdf", width=8, height=8)
    write.csv(distance_median,'between_group_distance_median.csv', row.names = T, quote = F)
  }
}

# correlation between taxa and pH
{
  keep = !is.na(metadata_all$sample_ph)
  reads_table = reads_table_all[,keep]
  metadata = metadata_all[keep,]
  
  reads_table = as.data.frame(t(reads_table))
  reads_table = reads_table + 0.5
  reads_table <- clr(reads_table)
  reads_table = as.data.frame(t(reads_table))
  
  line_list <- matrix(data = NA, nrow = nrow(reads_table), ncol =5)
  line_list <- as.data.frame(line_list)
  colnames(line_list) = c('Species','lm_p-value','lm_R-value','lm_formula','Slope_of_linear_regression')
  line_list$Species = row.names(reads_table)
  
  reads_table = reads_table[order(row.names(reads_table)),]

  p <- list()
  for (a in seq(nrow(reads_table))) {
    data = data.frame(pH = metadata$sample_ph, relative_abundance = as.numeric(as.character(reads_table[a,])))
    
    # linear
    linearMod = lm(pH~relative_abundance, data)
    
    pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])  # get p-value
    line_list[a,2]= pvalue
    r_squared <- summary(linearMod)$r.squared  # get R-value
    line_list[a,3]= r_squared
    
    b = format(round(linearMod$coefficients[2], 3), nsmall = 3) 
    line_list[a,5]= b
    c = format(round(linearMod$coefficients[1], 3), nsmall = 3)  
    if (linearMod$coefficients[1] > 0){  # output the formula of linear regression
      line_list[a,4]= paste0('y=',b,'x+',c)
    } else if (c < 0) {
      line_list[a,4]= paste0('y=',b,'x',c)
    } else {
      line_list[a,4]= paste0('y=',b,'x')
    }
    
    p[[a]] = ggplot(data, aes(x = pH, y = relative_abundance) ) +
        geom_point() +
        geom_smooth(method = "lm")+
        xlab("Vaginal pH")+
        ylab("Relative abundance")+
        labs(title = paste0(str_replace_all(row.names(reads_table)[a],"_"," "),'\n',line_list[a,4],'\n',"p-value = ",formatC(pvalue, format = "e", digits = 2),
                            '\n',"R-value = ",round(r_squared, 3)))+
        theme(axis.title = element_text(size = 12), 
            axis.text = element_text(size = 12), 
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 12))
#    ggsave(paste0("Fitting_",line_list$Species[a],'_model.pdf'), width=3, height=3)
    
  }
  
  library(ggpubr)
  ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],
            p[[11]],p[[12]],p[[13]],p[[14]],p[[15]],p[[16]],p[[17]],p[[18]],p[[19]],p[[20]],
            p[[21]],p[[22]],p[[23]],p[[24]],p[[25]],p[[26]],p[[27]],p[[28]],p[[29]],p[[30]],
            p[[31]],p[[32]],p[[33]],p[[34]],p[[35]],p[[36]],p[[37]],p[[38]],p[[39]],p[[40]])
  ggsave("Fitting_ph_linear_model.pdf", width=18.5, height=18.5)
  
  write.csv(line_list,'model_pH.csv')
  
}



########## impact of pH on the VMB #######
keep = !is.na(metadata_all$sample_ph)
reads_table_ph = reads_table_all[,keep]
metadata_ph = metadata_all[keep,]
mytypes_ph = mytypes[keep]

pic <- tsne$Y
pic = pic[keep,]

# classify pH
{
  x = metadata_ph$sample_ph
  x= (sort(x))
  data = data.frame(pH = x, sample = seq(length(x)))
  ggplot(data, aes(x=sample, y=pH)) + geom_point(size =0.2)+theme_bw()+ 
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  ggsave("ph_distribution.pdf", width=2, height=2)
  
  ph_0.33 = x[as.integer(length(x)/3)]
  ph_0.66 = x[as.integer(length(x)*2/3)]
  
  ph_0.33
  ph_0.66
  
  metadata_ph$sample_ph_new <- metadata_ph$sample_ph
  metadata_ph$sample_ph_new[metadata_ph$sample_ph <=ph_0.33] = 'Low'
  metadata_ph$sample_ph_new[metadata_ph$sample_ph >ph_0.66] = 'High'
  metadata_ph$sample_ph_new[metadata_ph$sample_ph >ph_0.33 & metadata_ph$sample_ph <=ph_0.66] = 'Medium'
  
}

# t-SNE
{
  # rarefy to normalize data
  reads_table_2 = as.data.frame(t(reads_table_ph))
  reads_table_2 = Rarefy(reads_table_2, depth = min(rowSums(reads_table_2)))
  reads_table_2 <- reads_table_2$otu.tab.rff
  reads_table_2 <- as.data.frame(reads_table_2)
  metadata_2 = as.data.frame(metadata_ph$sample_ph_new)

  pvalue = adonis2(reads_table_2 ~ ., data = metadata_2, method = "bray")
  pvalue
  
  
  reads_table = reads_table_ph
  reads_table = as.data.frame(t(reads_table))
  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
  reads_table <- reads_table$otu.tab.rff
  reads_table <- as.data.frame(reads_table)
  
  correlation <- adonis2(reads_table ~ sample_ph_new, data = metadata_ph, method = "bray")  # KitType significant
  output_pvalue <- "adonis_pvalue_pH" = correlation[1,5]
  
  pic <- data.frame(pic,metadata_ph$sample_ph_new)
  colnames(pic) <- c('X1','X2','pH')
  
  ggplot(pic, aes(X1, X2,color = pH))  +
    geom_point(size=0.2,aes(color = pH)) +
    xlab("t-SNE1") +
    ylab("t-SNE2") + 
    stat_ellipse(type = "t") + theme_bw()+
    theme(
      axis.title.x = element_text( size=7),
      axis.title.y = element_text( size=7),
      legend.text = element_text(size=7),
      legend.title = element_text(size=7)
    )
  ggsave("t_SNE_pH.pdf", width=3, height=2)
}

# abundance
{
  for (d in 1:3) {
    if (d==1) {
      keep = metadata_ph$sample_ph_new == 'Low' | metadata_ph$sample_ph_new == 'Medium'
    } else if (d==2) {
      keep = metadata_ph$sample_ph_new == 'Medium' | metadata_ph$sample_ph_new == 'High'
    } else {
      keep = metadata_ph$sample_ph_new == 'Low' | metadata_ph$sample_ph_new == 'High'
    }
    
    reads_table <- reads_table_ph[,keep]
    metadata <- as.character(metadata_ph$sample_ph_new[keep])
    metadata_factor = unique(as.character(metadata))
    
    if (d ==1) {
      p = dif_abundance(reads_table,metadata,paired_test = F,order_reverse =T, order = c('Medium',"Low"))
      p1 = as.data.frame(p$data[,c(1,5,13)])
      p1$Type = 'Low2Medium'
      p1$taxa = row.names(p1)
    } else if (d ==2) {
      p = dif_abundance(reads_table,metadata,paired_test = F, order = c('High',"Medium"))
      p2 = as.data.frame(p$data[,c(1,5,13)])
      p2$Type = 'Medium2High'
      p2$taxa = row.names(p2)
    } else {
      p = dif_abundance(reads_table,metadata,paired_test = F, order = c('High',"Low"))
      p3 = as.data.frame(p$data[,c(1,5,13)])
      p3$Type = 'Low2High'
      p3$taxa = row.names(p3)
    }
  }
  
  s = rbind(p1,p2,p3)
  
  colnames(s)[1] = 'Median difference in clr values'
  s = s[s$wi.eBH <= 0.05,]
  s = s[abs(s$`Median difference in clr values`) >= 1,]
  
  taxa_list = unique(s$taxa)
  taxa_order = data.frame("Taxa" = taxa_list,"Low2Medium"=NA, "Medium2High" = NA,"Low2High" = NA)
  for (j in 1:nrow(s)) {
    taxa_order[which(taxa_order$Taxa == s$taxa[j]), which(colnames(taxa_order) == s$Type[j])] <- s$`Median difference in clr values`[j]
  }
  taxa_order[is.na(taxa_order)] = 0
  taxa_order$delta = taxa_order$Low2Medium-taxa_order$Medium2High
  taxa_order = taxa_order[order(taxa_order$delta),]
  
  taxa_order$order2 = 0
  taxa_order$order2[taxa_order$Low2Medium > 0 ] = 1
  taxa_order$order2[taxa_order$Low2Medium < 0 ] = -1
  taxa_order = taxa_order[order(taxa_order$order2),]
  
  taxa_order$order3 = 0
  taxa_order$order3[taxa_order$Low2High > 0 ] = 1
  taxa_order$order3[taxa_order$Low2High < 0 ] = -1
  taxa_order = taxa_order[order(taxa_order$order3),]
  

  taxa_order$Taxa = str_replace_all(taxa_order$Taxa,'_'," ")
  s$taxa = str_replace_all(s$taxa,'_'," ")
  s$taxa <- factor(s$taxa, levels = taxa_order$Taxa)
  
  
  s$Type[s$Type == 'Low2Medium'] = "Low vs. Medium"
  s$Type[s$Type == 'Medium2High'] = "Medium vs. High"
  s$Type[s$Type == 'Low2High'] = "Low vs. High"
  
  s$Type <- factor(s$Type, levels = c( "Low vs. Medium","Medium vs. High","Low vs. High"))
  s$`Median difference in clr values`[s$`Median difference in clr values`>=5] =5
  
  ggplot(s, aes(taxa, Type)) + 
    geom_point(aes(col=`Median difference in clr values`, size=`-Log10(adj-pvalue)`)) + 
    scale_color_gradient2(midpoint=0, low="blue", mid="white",
                          high="red", space ="Lab" ) +
    coord_flip() +theme_bw()+          # convert x y axis
    labs(x = 'Taxa')+ 
    theme(axis.title = element_text(size = 12), 
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12))
  ggsave(paste('species_dif_abundance_ph.pdf',sep='_'),width=5.5, height=6.5)
}

# comparison of Lactobacillus vagitypes
{
  case_number = vector()
  case_number = c(sum(metadata_ph$Vagitype == 'Lactobacillus_crispatus'),
                  sum(metadata_ph$Vagitype == 'Lactobacillus_gasseri'),
                  sum(metadata_ph$Vagitype == 'Lactobacillus_jensenii'),
                  sum(metadata_ph$Vagitype == 'Lactobacillus_iners'),
                  sum(metadata_ph$sample_ph_new == 'Low'),
                  sum(metadata_ph$sample_ph_new == 'Medium'),
                  sum(metadata_ph$sample_ph_new == 'High'))
  case_number
  write.csv(case_number,'case_number_c_g_j_i_ph_L_M_H.csv')
  keep = metadata_ph$Vagitype == 'Lactobacillus_jensenii' | metadata_ph$Vagitype == 'Lactobacillus_crispatus'
  reads_table = reads_table_ph[,keep]
  metadata = metadata_ph[keep,]
  
  p <- dif_abundance(reads_table,metadata$Vagitype,order_reverse =T)    
  q5 = as.data.frame(p$data[,c(1,5,13)])
  q5$Type = 'crispatus_jensenii'
  q5$taxa = row.names(q5)
  
  keep = metadata_ph$Vagitype == 'Lactobacillus_iners' | metadata_ph$Vagitype == 'Lactobacillus_crispatus'
  reads_table = reads_table_ph[,keep]
  metadata = metadata_ph[keep,]
  
  p <- dif_abundance(reads_table,metadata$Vagitype,order_reverse =T)    
  q6 = as.data.frame(p$data[,c(1,5,13)])
  q6$Type = 'crispatus_iners'
  q6$taxa = row.names(q6)
  
  keep = metadata_ph$Vagitype == 'Lactobacillus_gasseri' | metadata_ph$Vagitype == 'Lactobacillus_crispatus'
  reads_table = reads_table_ph[,keep]
  metadata = metadata_ph[keep,]
  
  p <- dif_abundance(reads_table,metadata$Vagitype,order_reverse =T)    
  q7 = as.data.frame(p$data[,c(1,5,13)])
  q7$Type = 'crispatus_gasseri'
  q7$taxa = row.names(q7)
  
  
  keep = metadata_ph$Vagitype == 'Lactobacillus_jensenii' | metadata_ph$Vagitype == 'Lactobacillus_gasseri'
  reads_table = reads_table_ph[,keep]
  metadata = metadata_ph[keep,]
  
  p <- dif_abundance(reads_table,metadata$Vagitype,order_reverse =T)    
  q8 = as.data.frame(p$data[,c(1,5,13)])
  q8$Type = 'gasseri_jensenii'
  q8$taxa = row.names(q8)
  
  keep = metadata_ph$Vagitype == 'Lactobacillus_jensenii' | metadata_ph$Vagitype == 'Lactobacillus_iners'
  reads_table = reads_table_ph[,keep]
  metadata = metadata_ph[keep,]
  
  p <- dif_abundance(reads_table,metadata$Vagitype,order_reverse =F)    
  q9 = as.data.frame(p$data[,c(1,5,13)])
  q9$Type = 'jensenii_iners'
  q9$taxa = row.names(q9)
  
  keep = metadata_ph$Vagitype == 'Lactobacillus_iners' | metadata_ph$Vagitype == 'Lactobacillus_gasseri'
  reads_table = reads_table_ph[,keep]
  metadata = metadata_ph[keep,]
  
  p <- dif_abundance(reads_table,metadata$Vagitype,order_reverse =T)    
  q10 = as.data.frame(p$data[,c(1,5,13)])
  q10$Type = 'gasseri_iners'
  q10$taxa = row.names(q10)
  
  s = rbind(q5,q7,q6,q8,q9,q10)
  
  colnames(s)[1] = 'Median difference in clr values'
  s$Type[s$Type == "crispatus_iners"] = "L. crispatus vs. L. iners"
  s$Type[s$Type == "crispatus_gasseri"] = "L. crispatus vs. L. gasseri"
  s$Type[s$Type == "crispatus_jensenii"] = "L. crispatus vs. L. jensenii"
  s$Type[s$Type == "gasseri_iners"] = "L. gasseri vs. L. iners"
  s$Type[s$Type == "jensenii_iners"] = "L. jensenii vs. L. iners"
  s$Type[s$Type == "gasseri_jensenii"] = "L. gasseri vs. L. jensenii"
  
  s$Type = factor(s$Type, levels=c("L. crispatus vs. L. iners", "L. crispatus vs. L. gasseri", 
                                   "L. crispatus vs. L. jensenii","L. gasseri vs. L. iners",
                                   "L. jensenii vs. L. iners",'L. gasseri vs. L. jensenii'))
  
  colnames(s)[1] = 'Median difference in clr values'
  s = s[s$wi.eBH <= 0.05,]
  s = s[abs(s$`Median difference in clr values`) >= 1,]
  
  
  s$`Median difference in clr values`[s$`Median difference in clr values`>= 5] = 5
  s$`Median difference in clr values`[s$`Median difference in clr values`<= -5] = -5
  s$taxa = str_replace_all(s$taxa, "_", " ")
  s$taxa <- factor(s$taxa, levels = c('Ureaplasma',"Aerococcus","Gardnerella vaginalis","Atopobium vaginae",
                                      'Streptococcus',"Other Lactobacillus", "Lactobacillus iners",
                                      "Lactobacillus jensenii","Lactobacillus gasseri", "Lactobacillus crispatus"))
  
  ggplot(s, aes(taxa, Type)) + 
    geom_point(aes(col=`Median difference in clr values`, size=`-Log10(adj-pvalue)`)) + 
    scale_color_gradient2(midpoint=0, low="blue", mid="white",
                          high="red", space ="Lab" ) +
    coord_flip() +          # convert x y axis
    labs(x = 'Taxa')+ 
    theme(axis.title = element_text(size = 12), 
          axis.text = element_text(size = 12),
          axis.text.x = element_text(angle = 35, vjust = 1, hjust=1),
          legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12))  
  ggsave(paste('species_dif_abundance_within_lac.pdf',sep='_'),width=7, height=4)
  
}

############################## network  ############### 
Vagitype = vagitype(reads_table_all, th = 0.3)
metadata_all$Vagitype = Vagitype

# network all
{
  reads_table = reads_table_all
  row.names(reads_table) = str_replace_all(row.names(reads_table) , "_"," ")
  
  reads_table = as.data.frame(t(reads_table))
  reads_table = reads_table + 0.5
  reads_table <- clr(reads_table)      ### CLR normalization in rows
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
  
  paletteLength <- 50
  pvalue = 0.05;cor_parameter= 0;style = 1; pheatmap_fontsize = 5
  treeheight = 50; alpha = 0.05; mc.cores =60
  bar_max = max(p.yes.rr,na.rm = T)
  bar_min = min(p.yes.rr,na.rm = T)
  
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
  
  library(pheatmap)
  cluster = pheatmap(p.yes.rr, method = "complete",
                     clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean")
  
  taxa_cluster = sort(cutree(cluster$tree_row, k=6))
  taxa_cluster = as.data.frame(taxa_cluster)
  taxa_cluster$taxa = row.names(taxa_cluster)
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 1] = 11
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 4] = 12
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 5] = 13
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 6] = 14
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 2] = 15
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 3] = 16
  taxa_cluster$taxa_cluster = taxa_cluster$taxa_cluster-10
  cluster_color = cluster$tree_row$labels[cluster$tree_row$order]
  taxa_cluster = taxa_cluster[match(cluster_color,taxa_cluster$taxa),]
  write.csv(taxa_cluster,'network_all_taxa_cluster.csv')  
  taxa_cluster$taxa = NULL
  taxa_cluster$taxa_cluster = as.factor(taxa_cluster$taxa_cluster)
  colnames(taxa_cluster)[1] = "Taxa cluster"
  
  mycolors = list(`Taxa cluster` = c("1"="#FFDC00","2"="#aec7e8","3"="#17becf","4"="#e377c2","5"="#ff7f0e","6"="#d62728"))
  
  cluster = pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = pheatmap_fontsize, 
                     treeheight_row = treeheight, treeheight_col = treeheight, method = "complete",
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",show_colnames = T, 
                     annotation_colors = mycolors,annotation_row = taxa_cluster)
  
  pdf("network_all.pdf", width=6, height=5)
  print(cluster)
  dev.off()
  
}


## network for non-lac ph High

{
  keep = !is.na(metadata_all$sample_ph)
  reads_table = reads_table_all[,keep]
  metadata_ph = metadata_all[keep,]
  
  {
    x = metadata_ph$sample_ph
    x= (sort(x))
    data = data.frame(pH = x, sample = seq(length(x)))
    ggplot(data, aes(x=sample, y=pH)) + geom_point(size =0.2)+theme_bw()
    ggsave("ph_distribution.pdf", width=4, height=4)
    
    ph_0.33 = x[as.integer(length(x)/3)]
    ph_0.66 = x[as.integer(length(x)*2/3)]
    
    ph_0.33
    ph_0.66
    
    metadata_ph$sample_ph_new <- metadata_ph$sample_ph
    metadata_ph$sample_ph_new[metadata_ph$sample_ph <=ph_0.33] = 'Low'
    metadata_ph$sample_ph_new[metadata_ph$sample_ph >ph_0.66] = 'High'
    metadata_ph$sample_ph_new[metadata_ph$sample_ph >ph_0.33 & metadata_ph$sample_ph <=ph_0.66] = 'Medium'
    
  }
  
  keep = str_detect(metadata_ph$Vagitype,'Lactobacillus')
  reads_table_all_2 = reads_table[,!keep]
  metadata_all_2 = metadata_ph[!keep,]
  
  keep = metadata_all_2$sample_ph_new == 'High'
  cts_high <- reads_table_all_2[ ,keep]
  cts_high <- as.data.frame((cts_high))
  
  reads_table = cts_high
  row.names(reads_table) = str_replace_all(row.names(reads_table) , "_"," ")
  
  reads_table = as.data.frame(t(reads_table))
  reads_table = reads_table + 0.5
  reads_table <- clr(reads_table)      ### CLR normalization in rows
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
  
  paletteLength <- 50
  pvalue = 0.05;cor_parameter= 0;style = 1; pheatmap_fontsize = 5
  treeheight = 50; alpha = 0.05; mc.cores =60
  bar_max = max(p.yes.rr,na.rm = T)
  bar_min = min(p.yes.rr,na.rm = T)
  
  myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                seq(max(p.yes.rr,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))
  
  library(pheatmap)
  cluster = pheatmap(p.yes.rr, method = "complete",
                     clustering_distance_rows = "euclidean",clustering_distance_cols = "euclidean")
  
  taxa_cluster = sort(cutree(cluster$tree_row, k=6))
  taxa_cluster = as.data.frame(taxa_cluster)
  taxa_cluster$taxa = row.names(taxa_cluster)
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 1] = 11
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 4] = 12
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 6] = 13
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 6] = 14
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 2] = 13
  taxa_cluster$taxa_cluster[taxa_cluster$taxa_cluster == 3] = 15
  taxa_cluster$taxa_cluster = taxa_cluster$taxa_cluster-10
  cluster_color = cluster$tree_row$labels[cluster$tree_row$order]
  taxa_cluster = taxa_cluster[match(cluster_color,taxa_cluster$taxa),]
  write.csv(taxa_cluster,'network_all_taxa_cluster_high_non_lac.csv')  
  taxa_cluster$taxa = NULL
  taxa_cluster$taxa_cluster = as.factor(taxa_cluster$taxa_cluster)
  colnames(taxa_cluster)[1] = "Taxa cluster"
  
  mycolors = list(`Taxa cluster` = c("1"="#FFDC00","2"="#aec7e8","3"="#17becf","4"="#e377c2","5"="#ff7f0e","6"="#d62728"))
  
  cluster = pheatmap(p.yes.rr, color=myColor, breaks=myBreaks, fontsize = 6, 
                     treeheight_row = 20, treeheight_col = 20, method = "complete",
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",show_colnames = T)
  
  pdf("network_high_non_lac.pdf", width=5, height=5)
  print(cluster)
  dev.off()
  
}



########## BVT classification ############
setwd('/Users/binzhu/Desktop/pH')

data = read.csv('model_pH_network_MTS.csv',row.names = 1)
data_2 = data

data[data=='a'] = 1
data[data=='b'] = 2
data[data=='c'] = 3
data[data=='d'] = 4
data[data=='e'] = 5
data[data=='f'] = 6
data[data=='g'] = 7
data[data=='h'] = 8
data[data=='i'] = 9
data <- data.frame(apply(data, 2, function(x) as.numeric(as.character(x))))
row.names(data) = row.names(data_2)

# https://medium.com/@maryam.alizadeh/clustering-categorical-or-mixed-data-in-r-c0fb6ff38859
library(cluster)
df<-data_2
df$Sensitivity_to_pH = as.factor(df$Sensitivity_to_pH)
df$Taxa_correlation_network = as.factor(df$Taxa_correlation_network)
df$Interaction_with_host = as.factor(df$Interaction_with_host)
d_dist<-daisy(df, metric = "gower")

hcl <- hclust(d_dist, method = "complete")
output = hcl$labels[hcl$order]
#write.csv(output,'Cluster.csv')  # add cluster group mannually 

hcl$labels = str_replace_all(hcl$labels,'_',' ')

library(shipunov)
pdf("Taxa_classification_2.pdf", width=12, height=6)
old.par <- par(mfrow=c(1, 2))
plot(hcl, labels=gsub("[A-z.]", "  ", hcl$labels))
Tctext(hcl, srt=90, add=0.04, adj=c(1, 0.5))
plot(hcl)
par(old.par)
dev.off()

pdf("Taxa_classification_1.pdf", width=12, height=6)
old.par <- par(mfrow=c(1, 2))
plot(hcl, labels=gsub("[A-z.]", "  ", hcl$labels))
Tctext(hcl, srt=90, add=0.04, adj=c(1, 0.5))
plot(hcl)
rect.hclust(hcl, k=7, border="red")
par(old.par)
dev.off()


library(pheatmap)
taxa_order = read.csv('Cluster.csv')
taxa_order = taxa_order[order(taxa_order$Taxa),]
taxa_order = taxa_order[order(taxa_order$Cluster),]
row.names(taxa_order) = taxa_order$Taxa
taxa_order$Taxa = NULL
row.names(taxa_order) = str_replace_all(row.names(taxa_order), '_'," ")
row.names(data) = str_replace_all(row.names(data), '_'," ")


data_2 = data[match(row.names(taxa_order),row.names(data)),]
colnames(data_2) =c("Host response", "Taxa correlation", "Vaginal pH")
data_2 = data_2[c( "Taxa correlation", "Vaginal pH","Host response")]

#mycolors = list(Cluster = c("I"="#FFDC00","II"="#49C800","III"="#aec7e8","IV"="#17becf","V"="#7f7f7f",
#                            "VI"="#ff7f0e","VII"="#FF00D1","VIII"="#d62728","Others"="#9C9146"))
mycolors = list(Cluster = c("I"="#FFDC00","II"="#49C800","III"="#aec7e8","IV"="#17becf","V"="#7f7f7f",
                            "VI"="#ff7f0e","VII"="#d62728","Others"="#9C9146"))

cluster = pheatmap(data_2, cluster_rows=F, show_rownames=TRUE,
                   cluster_cols=F, annotation_row = taxa_order,annotation_colors =mycolors)

pdf("Taxa_classification_3.pdf", width=3.3, height=5)
print(cluster)
dev.off()


########## application of B-vagitype ############
metadata = metadata_all[,str_detect(colnames(metadata_all), 
                                    'dx_bv|sample_ph|african_american|caucasian|education|vodor|vdischarge')]


# get vagitype
{
  reads_table_all_abundance = get_abundance_table(reads_table_all)
  type_th = 0.3
  
  reads_table = as.data.frame(t(reads_table_all_abundance))
  reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
  mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
  maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
  mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
  mytypes[maxprop < type_th] <- "No_Type"
  uniqueTypes <- unique(mytypes)
  
  top_vagitype = as.data.frame(table(mytypes))
  top_vagitype = top_vagitype[order(top_vagitype$Freq, decreasing = T),]
}

# get biological vagitype (B-vagitype)
{
  BVT = as.data.frame(mytypes)
  taxa_order = read.csv('Cluster.csv')
#  taxa_order[nrow(taxa_order)+1,2] = 'IX'
#  taxa_order[nrow(taxa_order),1] = 'No_Type'
  
  n = match(BVT$mytypes, taxa_order$Taxa)
  BVT$BVT = taxa_order$Cluster[n]
  BVT = BVT$BVT
  BVT[is.na(BVT)] = 'Others'
  keep = mytypes %in% top_vagitype$mytypes[1:12]
  mytypes[!keep] = 'Others'
}


# get CST (community state types) 1
{
  reads_table_cluster <- as.matrix(reads_table_all_abundance)
  Dominant_species_2 <- as.data.frame(BVT)
  Dominant_species_2$Vagitype <- mytypes
  row.names(Dominant_species_2) = colnames(reads_table_cluster)
  
#  reads_table = as.data.frame(t(reads_table_all))
#  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
#  reads_table <- reads_table$otu.tab.rff
#  reads_table <- as.data.frame(reads_table)
#  Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
#  Bray_Curtis = as.dist(Bray_Curtis)
  
  CST = pheatmap(reads_table_cluster,cluster_rows = T, cluster_cols =T, 
                 clustering_distance_cols = "euclidean", clustering_method = "complete")
  
  CST_cluster = sort(cutree(CST$tree_col, k=9))
  CST_cluster = as.data.frame(CST_cluster)
  
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 1] = 'I'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 2] = 'III'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 3] = 'II'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 4] = 'IVC'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 5] = 'IVA'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 6] = 'IVB'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 7] = 'IVD'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 8] = 'V'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 9] = 'IVE'

  Dominant_species_2$CST = CST_cluster$CST_cluster[match(row.names(Dominant_species_2), row.names(CST_cluster))]
  Dominant_species_2$CST = as.factor(Dominant_species_2$CST)
  
  
  row.names(reads_table_cluster) = str_replace_all(row.names(reads_table_cluster),'_',' ')
  Dominant_species_2$Vagitype = str_replace_all(Dominant_species_2$Vagitype ,'_',' ')
  Dominant_species_2$Vagitype = as.factor(Dominant_species_2$Vagitype)
  write.csv(Dominant_species_2, 'VMB_type.csv')
  
  y = table(Dominant_species_2$BVT)
  write.csv(y, 'VMB_type_number_BVT.csv')
  
  y = table(Dominant_species_2$CST)
  write.csv(y, 'VMB_type_number_CST.csv')
  
  mycolors = list(Vagitype = c("Lactobacillus crispatus"="#FFDC00","Lactobacillus gasseri"="#c5b0d5","Lactobacillus iners"="#aec7e8",
                               "Lactobacillus jensenii"="#49C800","Ca. L. vaginae"="#ff7f0e","Gardnerella vaginalis"="#d62728","No Type"="#4d4d4d",
                               "Atopobium vaginae"="#c49c94","Prevotella bivia"="#17becf","Mycoplasma"="#e377c2",
                               "Sneathia amnii"="#FF00D1", "Streptococcus"="#7f7f7f","Others"="#9C9146"),
                  CST = c("I"="#FFDC00","II"="#c5b0d5","III"="#aec7e8","V"="#49C800","IVA"="#ff7f0e","IVB"="#d62728","IVC"="#4d4d4d","IVD"="#17becf",
                               "IVE"="#e377c2"),
                  BVT = c("I"="#FFDC00","II"="#49C800","III"="#aec7e8","VI"="#ff7f0e","VII"="#d62728","IV"="#17becf","V"="#7f7f7f",
                          "Others"="#9C9146"))
  
  pdf("heatmap_clustering.pdf", width=18, height=8)
  print(pheatmap(reads_table_cluster,cluster_rows = T, cluster_cols =T, clustering_distance_cols = "euclidean", clustering_method = "complete",
                 annotation_col = Dominant_species_2, annotation_colors =mycolors,show_colnames = F,treeheight_col = 50))
  dev.off()
  
  CST_cluster$SampleID = row.names(CST_cluster)
  CST_cluster = CST_cluster[match(row.names(metadata_all),CST_cluster$SampleID),]
  CST = CST_cluster$CST_cluster
  
#  keep = mytypes %in% top_vagitype$mytypes[1:8]
#  mytypes[!keep] = 'Others'
  Dominant_species_2$Vagitype = mytypes
}

metadata_test = data.frame(CST = CST, Vagitype= mytypes, BVT = BVT) 

for (n in 1:3) {
  output_table = as.data.frame(matrix(data = NA, ncol = 7, nrow = length(unique(metadata_test[,n]))))

  colnames(output_table) = colnames(metadata)
  row.names(output_table) = sort(unique(metadata_test[,n]))
  
  output_table_2 = as.data.frame(matrix(data = NA, ncol = 7, nrow = 2))
  colnames(output_table_2) = colnames(metadata)
  row.names(output_table_2) = c('Kruskal–Wallis test',"Pearson's Chi-squared test")
  
  
  for (y in 1: ncol(output_table)) {
    for (x in 1:nrow(output_table)) {
      if (y == 3 | y == 6) {
        data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
        data = data[!is.na(data)]
        data_2 = mean(data)
        data_2 = round(data_2, 1)
        data_3 = sd(data)
        data_3 = round(data_3, 1)
        output_table[x,y] = paste0(data_2, " (",data_3, "); ",length(data))
      } else {
        data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
        data = data[!is.na(data)]
        data_2 = sum(data == 'Yes') / length(data)
        data_2 = round(data_2, 3)
        output_table[x,y] = paste0(data_2*100,'%; ',length(data))
      }
    }
    data_4 = as.data.frame(cbind(metadata_test[,n], metadata[,y]))
    keep = !is.na(data_4[,2])
    data_4 = data_4[keep,]
    output_table_2[1,y] = kruskal.test(V2~V1, data_4)$p.value
    output_table_2[2,y] = chisq.test(data_4$V2,  data_4$V1, correct=FALSE)$p.value
  }
  output_table = rbind(output_table,output_table_2)
  write.csv(output_table, paste0('Distribution_',colnames(metadata_test)[n],'.csv'))
}





########## application of B-vagitype 5 types ############
metadata = metadata_all[,str_detect(colnames(metadata_all), 
                                    'dx_bv|sample_ph|african_american|caucasian|education|vodor|vdischarge')]


# get vagitype
{
  reads_table_all_abundance = get_abundance_table(reads_table_all)
  type_th = 0.3
  
  reads_table = as.data.frame(t(reads_table_all_abundance))
  reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
  mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
  maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
  mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
  mytypes[maxprop < type_th] <- "No_Type"
  uniqueTypes <- unique(mytypes)
  
  top_vagitype = as.data.frame(table(mytypes))
  top_vagitype = top_vagitype[order(top_vagitype$Freq, decreasing = T),]
}

# get biological vagitype (B-vagitype)
{
  BVT = as.data.frame(mytypes)
  taxa_order = read.csv('Cluster.csv')

  n = match(BVT$mytypes, taxa_order$Taxa)
  BVT$BVT = taxa_order$Cluster[n]
  BVT = BVT$BVT
  BVT[is.na(BVT)] = 'Others'
  BVT[BVT != "I" & BVT != "VI" & BVT != "III" & BVT != "VII" ] = 'Others'
}

keep = mytypes %in% top_vagitype$mytypes[1:4]
mytypes[!keep] = 'Others'

# get CST (community state types) 1
{
  reads_table_cluster <- as.matrix(reads_table_all_abundance)
  Dominant_species_2 <- as.data.frame(BVT)
  Dominant_species_2$Vagitype <- mytypes
  row.names(Dominant_species_2) = colnames(reads_table_cluster)
  
#  reads_table = as.data.frame(t(reads_table_all))
#  reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
#  reads_table <- reads_table$otu.tab.rff
#  reads_table <- as.data.frame(reads_table)
#  Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
#  Bray_Curtis = as.dist(Bray_Curtis)
  
  CST = pheatmap(reads_table_cluster,cluster_rows = T, cluster_cols =T, 
                 clustering_distance_cols = "euclidean", clustering_method = "complete")
  
  CST_cluster = sort(cutree(CST$tree_col, k=10))
  CST_cluster = as.data.frame(CST_cluster)

  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 1] = 'I'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 2] = 'III'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 3] = 'II'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 4] = 'IVC'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 5] = 'IVA'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 6] = 'IVB'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 7] = 'IVD'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 8] = 'V'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 9] = 'IVE'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 10] = 'IVD'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster != 'I' & CST_cluster$CST_cluster != 'III' &
                            CST_cluster$CST_cluster != 'IVA' &CST_cluster$CST_cluster != 'IVB' ] = "Others"
  
  Dominant_species_2$CST_2 = CST_cluster$CST_cluster[match(row.names(Dominant_species_2), row.names(CST_cluster))]
  Dominant_species_2$CST_2 = as.factor(Dominant_species_2$CST_2)
  
  
  
  CST_cluster = sort(cutree(CST$tree_col, k=5))
  CST_cluster = as.data.frame(CST_cluster)
  
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 1] = 'I'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 2] = 'III'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 3] = 'II'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 4] = 'IV'
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 5] = 'V'
 
  Dominant_species_2$CST_1 = CST_cluster$CST_cluster[match(row.names(Dominant_species_2), row.names(CST_cluster))]
  Dominant_species_2$CST_1 = as.factor(Dominant_species_2$CST_1)

  row.names(reads_table_cluster) = str_replace_all(row.names(reads_table_cluster),'_',' ')
  Dominant_species_2$Vagitype = str_replace_all(Dominant_species_2$Vagitype ,'_',' ')

  Dominant_species_2$Vagitype = factor(Dominant_species_2$Vagitype, levels = c("Lactobacillus crispatus",
                                                                               "Lactobacillus iners",
                                                                               "Gardnerella vaginalis",
                                                                               "Ca. L. vaginae", "Others"))
  
  Dominant_species_2$BVT = factor(Dominant_species_2$BVT)
  colnames(Dominant_species_2)[1] = 'BVT'
  
  
  mycolors = list(Vagitype = c("Lactobacillus crispatus"="#FFDC00","Lactobacillus iners"="#aec7e8","Ca. L. vaginae"="#ff7f0e",
                               "Gardnerella vaginalis"="#d62728","Others"="#9C9146"),
                  CST_2 = c("I"="#FFDC00","III"="#aec7e8","IVA"="#ff7f0e","IVB"="#d62728","Others"="#9C9146"),
                  CST_1 = c("I"="#FFDC00","II"="#c5b0d5","III"="#aec7e8","V"="#49C800","IV"="#4d4d4d"),
                  BVT = c("I"="#FFDC00","III"="#aec7e8","VI"="#ff7f0e","VII"="#d62728","Others"="#9C9146"))
  
  
  pdf("heatmap_clustering_5_CST_5_vagitype.pdf", width=18, height=8)
  print(pheatmap(reads_table_cluster,cluster_rows = T, cluster_cols =T, clustering_distance_cols = "euclidean", clustering_method = "complete",
                 annotation_col = Dominant_species_2,show_colnames = F,treeheight_col = 50, annotation_colors =mycolors))
  dev.off()
  
  
}

# compare
metadata_test = Dominant_species_2
for (n in 1:4) {
  output_table = as.data.frame(matrix(data = NA, ncol = 7, nrow = 5))
  colnames(output_table) = colnames(metadata)
  row.names(output_table) = sort(unique(metadata_test[,n]))
  
  output_table_2 = as.data.frame(matrix(data = NA, ncol = 7, nrow = 2))
  colnames(output_table_2) = colnames(metadata)
  row.names(output_table_2) = c('Kruskal_Wallis test',"Pearson's Chi-squared test")
  
  for (y in 1: ncol(output_table)) {
    for (x in 1:nrow(output_table)) {
      if (y == 3 | y == 6) {
        data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
        data = data[!is.na(data)]
        data_2 = median(data)
        data_2 = round(data_2, 1)
        data_3 = quantile(data)
        data_3 = round(data_3, 1)
        output_table[x,y] = paste0(data_2, " (",data_3[2], ", ",data_3[4], "); ",length(data))
      } else {
        data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
        data = data[!is.na(data)]
        data_2 = sum(data == 'Yes') / length(data)
        data_2 = round(data_2, 3)
        output_table[x,y] = paste0(data_2*100,'%; ',length(data))
      }
    }
    
    data_4 = as.data.frame(cbind(metadata_test[,n], metadata[,y]))
    keep = !is.na(data_4[,2])
    data_4 = data_4[keep,]
    output_table_2[1,y] = kruskal.test(V2~V1, data_4)$p.value
    output_table_2[2,y] = chisq.test(data_4$V2,  data_4$V1, correct=FALSE)$p.value
    
  }
  output_table = rbind(output_table , output_table_2)
  write.csv(output_table, paste0('Distribution_5_types_',colnames(metadata_test)[n],'.csv'))
}
for (n in 1:4) {
  output_table = as.data.frame(matrix(data = NA, ncol = 7, nrow = 5))
  colnames(output_table) = colnames(metadata)
  row.names(output_table) = sort(unique(metadata_test[,n]))
  
  output_table_2 = as.data.frame(matrix(data = NA, ncol = 7, nrow = 2))
  colnames(output_table_2) = colnames(metadata)
  row.names(output_table_2) = c('Kruskal_Wallis test',"Pearson's Chi-squared test")
  
  for (y in 1: ncol(output_table)) {
    for (x in 1:nrow(output_table)) {
      if (y == 3 | y == 6) {
        data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
        data = data[!is.na(data)]
        data_2 = mean(data)
        data_2 = round(data_2, 1)
        output_table[x,y] = data_2
      } else {
        data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
        data = data[!is.na(data)]
        data_2 = sum(data == 'Yes') / length(data)
        data_2 = round(data_2, 3)
        output_table[x,y] = data_2
      }
    }
    
    data_4 = as.data.frame(cbind(metadata_test[,n], metadata[,y]))
    keep = !is.na(data_4[,2])
    data_4 = data_4[keep,]
    output_table_2[1,y] = kruskal.test(V2~V1, data_4)$p.value
    output_table_2[2,y] = chisq.test(data_4$V2,  data_4$V1, correct=FALSE)$p.value
    
  }
  output_table = rbind(output_table , output_table_2)
  write.csv(output_table, paste0('Distribution_5_types_2_',colnames(metadata_test)[n],'.csv'))
}



########## lower sample size and application of B-vagitype ############
compare_results_1 = as.data.frame(matrix(data = NA, ncol = 8, nrow = 0))
colnames(compare_results_1) = c(colnames(metadata),'Classification')

compare_results_2 = as.data.frame(matrix(data = NA, ncol = 8, nrow = 0))
colnames(compare_results_2) = c(colnames(metadata),'Classification')

repeat_time = 100
sample_per_time = 100
options(show.error.locations = TRUE)
for (x  in 1: repeat_time) {
    set.seed(x)
    keep = sample(c(1:nrow(metadata_all)), sample_per_time)
    
    metadata_all_2 = metadata_all[keep,]
    reads_table_all_2 = reads_table_all[,keep]
    
    metadata = metadata_all_2[,str_detect(colnames(metadata_all_2), 
                                          'dx_bv|sample_ph|african_american|caucasian|education|vodor|vdischarge')]
    
    
    # get vagitype
    {
      reads_table_all_abundance = get_abundance_table(reads_table_all_2)
      type_th = 0.3
      
      reads_table = as.data.frame(t(reads_table_all_abundance))
      reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
      mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
      maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
      mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
      mytypes[maxprop < type_th] <- "No_Type"
      uniqueTypes <- unique(mytypes)
      
      top_vagitype = as.data.frame(table(mytypes))
      top_vagitype = top_vagitype[order(top_vagitype$Freq, decreasing = T),]
      top_vagitype
      keep = mytypes %in% top_vagitype$mytypes[1:12]
      mytypes[!keep] = 'Others'
      unique(mytypes)
    }
    
    # get CST (community state types)
    {
      reads_table_cluster <- as.matrix(reads_table_all_abundance)
      Dominant_species_2 <- as.data.frame(mytypes)
      colnames(Dominant_species_2) = 'Vagitype'
      row.names(Dominant_species_2) = colnames(reads_table_cluster)
      
      mycolors = list(Vagitype = c(Lactobacillus_iners="#aec7e8",Lactobacillus_crispatus="#FFDC00",
                                   Gardnerella_vaginalis="#d62728",Ca._L._vaginae="#ff7f0e",Others="#faf0e6"))
      
      
      # bray-curtis
      # rarefy to normalize data
      reads_table = as.data.frame(t(reads_table_all_2))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table <- reads_table$otu.tab.rff
      reads_table <- as.data.frame(reads_table)
      Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
      Bray_Curtis = as.dist(Bray_Curtis)
      
      CST = pheatmap(reads_table_cluster,cluster_rows = T, cluster_cols =T, clustering_distance_cols = Bray_Curtis, 
                     clustering_method = "ward.D")
      
      CST_cluster = sort(cutree(CST$tree_col, k=5))
      CST_cluster = as.data.frame(CST_cluster)
      #  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 1] = 'IVB'
      #  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 2] = 'I'
      #  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 3] = 'III'
      #  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 4] = 'IVA'
      #  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 5] = 'Others'
      
      Dominant_species_2$CST = CST_cluster$CST_cluster[match(row.names(Dominant_species_2), row.names(CST_cluster))]
      Dominant_species_2$CST = as.factor(Dominant_species_2$CST)
      
      CST_cluster$SampleID = row.names(CST_cluster)
      CST_cluster = CST_cluster[match(row.names(metadata_all_2),CST_cluster$SampleID),]
      CST = CST_cluster$CST_cluster
      CST_1 = CST
      unique(CST)
    }
    
    # get CST (community state types) 2
    {
      reads_table_cluster <- as.matrix(reads_table_all_abundance)
      Dominant_species_2 <- as.data.frame(mytypes)
      colnames(Dominant_species_2) = 'Vagitype'
      row.names(Dominant_species_2) = colnames(reads_table_cluster)
      
      mycolors = list(Vagitype = c(Lactobacillus_iners="#aec7e8",Lactobacillus_crispatus="#FFDC00",
                                   Gardnerella_vaginalis="#d62728",Ca._L._vaginae="#ff7f0e",Others="#faf0e6"))
      
      
      # bray-curtis
      # rarefy to normalize data
      reads_table = as.data.frame(t(reads_table_all_2))
      reads_table = Rarefy(reads_table, depth = min(rowSums(reads_table)))
      reads_table <- reads_table$otu.tab.rff
      reads_table <- as.data.frame(reads_table)
      Bray_Curtis <- as.matrix(vegdist(reads_table, method = "bray", binary=FALSE))
      Bray_Curtis = as.dist(Bray_Curtis)
      
      CST = pheatmap(reads_table_cluster,cluster_rows = T, cluster_cols =T, clustering_distance_cols = Bray_Curtis, 
                     clustering_method = "ward.D")
      
      CST_cluster = sort(cutree(CST$tree_col, k=10))
      CST_cluster = as.data.frame(CST_cluster)
      CST_cluster_result = as.data.frame(table(CST_cluster$CST_cluster))
      CST_cluster_result = CST_cluster_result[order(CST_cluster_result$Freq, decreasing = T),]
      CST_cluster_result = as.numeric(as.character(CST_cluster_result$Var1[1:4]))
      CST_cluster$CST_cluster[!CST_cluster$CST_cluster %in% CST_cluster_result] = 'Others'
      
      Dominant_species_2$CST = CST_cluster$CST_cluster[match(row.names(Dominant_species_2), row.names(CST_cluster))]
      Dominant_species_2$CST = as.factor(Dominant_species_2$CST)
      
      CST_cluster$SampleID = row.names(CST_cluster)
      CST_cluster = CST_cluster[match(row.names(metadata_all_2),CST_cluster$SampleID),]
      CST = CST_cluster$CST_cluster
      CST_2 = CST
      unique(CST)
    }
    
    # get biological vagitype (B-vagitype)
    {
      taxa_order = read.csv('Cluster.csv')
      taxa_order = taxa_order[order(taxa_order$Taxa),]
      taxa_order = taxa_order[order(taxa_order$Cluster),]
      row.names(taxa_order) = taxa_order$Taxa
      taxa_order$Taxa = NULL
      
      type_th = 0.3
      reads_table = as.data.frame(t(reads_table_all_abundance))
      reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
      mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
      maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
      mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
      mytypes[maxprop < type_th] <- "No_Type"
      
      BVT = as.data.frame(mytypes)
      taxa_order = as.data.frame(taxa_order)
      taxa_order$Taxa = row.names(taxa_order)
      
      n = match(BVT$mytypes, taxa_order$Taxa)
      BVT$BVT = taxa_order$Cluster[n]
      BVT = BVT$BVT
      table(BVT)
      BVT[is.na(BVT)] = 'Others'
      BVT[BVT != 'I' & BVT != 'VI' & BVT != 'III' & BVT != 'VII'] = 'Others'
    }
    
    #  comparison
    {
      # get vagitype
      {
        reads_table = as.data.frame(t(reads_table_all_abundance))
        reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
        mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
        maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
        mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
        mytypes[maxprop < type_th] <- "No_Type"
        uniqueTypes <- unique(mytypes)
        
        top_vagitype = as.data.frame(table(mytypes))
        top_vagitype = top_vagitype[order(top_vagitype$Freq, decreasing = T),]
        top_vagitype
        keep = mytypes %in% top_vagitype$mytypes[1:4]
        mytypes[!keep] = 'Others'
        unique(mytypes)
      }
      
      # compare
      metadata_test = data.frame(CST_1 = CST_1,CST_2 = CST_2, Vagitype= mytypes, BVT = BVT) 
      
      for (n in 1:4) {
        output_table = as.data.frame(matrix(data = NA, ncol = 7, nrow = 5))
        colnames(output_table) = colnames(metadata)
        row.names(output_table) = sort(unique(metadata_test[,n]))
        
        output_table_2 = as.data.frame(matrix(data = NA, ncol = 7, nrow = 2))
        colnames(output_table_2) = colnames(metadata)
        row.names(output_table_2) = c('Kruskal_Wallis test',"Pearson's Chi-squared test")
        
        for (y in 1: ncol(output_table)) {
          for (x in 1:nrow(output_table)) {
            if (y == 3 | y == 6) {
              data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
              data = data[!is.na(data)]
              data_2 = mean(data)
              data_2 = round(data_2, 1)
              data_3 = sd(data)
              data_3 = round(data_3, 1)
              output_table[x,y] = paste0(data_2, " (",data_3, "); ",length(data))
            } else {
              data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
              data = data[!is.na(data)]
              data_2 = sum(data == 'Yes') / length(data)
              data_2 = round(data_2, 3)
              output_table[x,y] = paste0(data_2*100,'%; ',length(data))
            }
          }
          
          data_4 = as.data.frame(cbind(metadata_test[,n], metadata[,y]))
          keep = !is.na(data_4[,2])
          data_4 = data_4[keep,]
          output_table_2[1,y] = kruskal.test(V2~V1, data_4)$p.value
          
          if (length(unique(data_4$V2)) < 2) {
            next
          }
          output_table_2[2,y] = chisq.test(data_4$V2,  data_4$V1, correct=FALSE)$p.value
          
          output_table_2$Classification = colnames(metadata_test)[n]
        }
        compare_results_1 = rbind(compare_results_1,output_table_2[1,])
        compare_results_2 = rbind(compare_results_2,output_table_2[2,])
      }
      
    }
    

  
  
}

Classification = compare_results_1$Classification
compare_results_1$Classification = NULL
compare_results_2$Classification = NULL

compare_results_1 = -log10(compare_results_1)
compare_results_2 = -log10(compare_results_2)

compare_results_1$Classification = Classification
compare_results_2$Classification = Classification

variables_list = colnames(metadata)

for (a in 1:length(variables_list)) {

    data = compare_results_2[,c(a,8)]
    colnames(data)[1] = 'Value'
    data$Classification = factor(data$Classification, levels = c("BVT","CST_1","CST_2","Vagitype"))

    # calculate significance
    factor_levels = c("BVT","CST_1","CST_2","Vagitype")
    n = length(factor_levels)
    
    sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
    colnames(sig) = factor_levels
    row.names(sig) = factor_levels

    for (c in 1:(n-1)) {
      for (b in (c+1) : n) {
        factor_level1 <- data[data$Classification == factor_levels[c],1]
        factor_level2 <- data[data$Classification == factor_levels[b],1]
        
        sig[c,b] <- wilcox.test(factor_level1, factor_level2, paired = T)$p.value
      }
    }
    write.csv(sig,paste0('100_random_',variables_list[a],'.csv'))
    
    data$Value[data$Value > 5] = 5
    ggplot(data=data, aes(x=Classification, y=Value)) +
      geom_boxplot(fill='gray', color="black", outlier.shape=NA) +
      geom_jitter(shape=16, size = 0.1)+theme_bw()+
      ylab(paste0(variables_list[a],"\nPearson's Chi-squared test\n-Log10(p-values)"))+
      theme(axis.title = element_text(size = 9), 
            axis.text = element_text(size = 9), 
            legend.text = element_text(size = 9), 
            legend.title = element_text(size = 9),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    ggsave(paste0('100_random_',variables_list[a],'.pdf'),width=1.7, height=3)
}





########## ####### pick up samples for MTS ###### 
setwd('/Users/binzhu/Desktop/pH/results')
reads_table_all_abundance = get_abundance_table(reads_table_all)
type_th = 0.3

# get vagitype
{
  reads_table = as.data.frame(t(reads_table_all_abundance))
  reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
  
  mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
  maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
  mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
  mytypes[maxprop < type_th] <- "No_Type"
  uniqueTypes <- unique(mytypes)
  
  data = cbind(colnames(reads_table_all), mytypes)
  write.csv(data,'Vagitype.csv')
  top_vagitype = as.data.frame(table(mytypes))
  top_vagitype = top_vagitype[order(top_vagitype$Freq, decreasing = T),]
  keep = mytypes %in% top_vagitype$mytypes[1:13]
  
  mytypes[!keep] = 'Others'
  metadata_all$Vagitype = mytypes
}

keep = mytypes != "Others"
metadata_all_2 = metadata_all[keep,]
reads_table_all_abundance_2 = reads_table_all_abundance[,keep]
reads_table_all_abundance_2 = as.data.frame(reads_table_all_abundance_2)
mytypes = mytypes[keep]
top_10 = unique(mytypes)

keep = metadata_all_2$dx_bv == 'No' & !is.na(metadata_all_2$dx_bv) & 
  metadata_all_2$dx_warts == 'No' & !is.na(metadata_all_2$dx_warts) & 
  metadata_all_2$dx_UTI == 'No' & !is.na(metadata_all_2$dx_UTI) & 
  metadata_all_2$dx_trich == 'No' & !is.na(metadata_all_2$dx_trich) & 
  metadata_all_2$dx_yeast_infection == 'No' & !is.na(metadata_all_2$dx_yeast_infection) & 
  metadata_all_2$yeast == 'No' & !is.na(metadata_all_2$yeast) & 
  metadata_all_2$bv == 'No' & !is.na(metadata_all_2$bv) & 
  metadata_all_2$vodor == 'No' & !is.na(metadata_all_2$vodor) & 
  metadata_all_2$vitching == 'No' & !is.na(metadata_all_2$vitching) & 
  metadata_all_2$vdischarge == 'No' & !is.na(metadata_all_2$vdischarge) & 
  metadata_all_2$quinolones_status == 'No' & !is.na(metadata_all_2$quinolones_status) & 
  metadata_all_2$sulfa_status == 'No' & !is.na(metadata_all_2$sulfa_status) & 
  metadata_all_2$penicillins_status == 'No' & !is.na(metadata_all_2$penicillins_status) & 
  metadata_all_2$other_antibiotics_status == 'No' & !is.na(metadata_all_2$other_antibiotics_status) & 
  metadata_all_2$cephalosporins_status == 'No' & !is.na(metadata_all_2$cephalosporins_status) & 
  metadata_all_2$macrolides_status == 'No' & !is.na(metadata_all_2$macrolides_status) & 
  metadata_all_2$doctor_yeast_infection_medication_status == 'No' & !is.na(metadata_all_2$doctor_yeast_infection_medication_status) & 
  metadata_all_2$nitrofurantoin_status == 'No' & !is.na(metadata_all_2$nitrofurantoin_status) & 
  metadata_all_2$insulin_status == 'No' & !is.na(metadata_all_2$insulin_status) & 
  metadata_all_2$oral_steroids_status == 'No' & !is.na(metadata_all_2$oral_steroids_status) & 
  metadata_all_2$tetracyclines_status == 'No' & !is.na(metadata_all_2$tetracyclines_status) & 
  metadata_all_2$counter_yeast_infection_medication_status == 'No' & !is.na(metadata_all_2$counter_yeast_infection_medication_status) & 
  metadata_all_2$urinary_tract_infection_status == 'No' & !is.na(metadata_all_2$urinary_tract_infection_status) & 
  metadata_all_2$metrodinazole_status == 'No' & !is.na(metadata_all_2$metrodinazole_status) & 
  metadata_all_2$syphilis_status == 'No' & !is.na(metadata_all_2$syphilis_status) & 
  metadata_all_2$gonorrhea_status == 'No' & !is.na(metadata_all_2$gonorrhea_status) & 
  metadata_all_2$yeast_infection_status == 'No' & !is.na(metadata_all_2$yeast_infection_status) & 
  metadata_all_2$genital_warts_status == 'No' & !is.na(metadata_all_2$genital_warts_status) & 
  metadata_all_2$trichomoniasis_status == 'No' & !is.na(metadata_all_2$trichomoniasis_status) & 
  metadata_all_2$hpv_status == 'No' & !is.na(metadata_all_2$hpv_status) & 
  metadata_all_2$pelvic_inflammatory_disease_status == 'No' & !is.na(metadata_all_2$pelvic_inflammatory_disease_status) & 
  metadata_all_2$genital_herpes_status == 'No' & !is.na(metadata_all_2$genital_herpes_status) & 
  metadata_all_2$chlamydia_status == 'No' & !is.na(metadata_all_2$chlamydia_status) & 
  metadata_all_2$abnormal_pap_smear_status == 'No' & !is.na(metadata_all_2$abnormal_pap_smear_status) & 
  metadata_all_2$vaginitis_status == 'No' & !is.na(metadata_all_2$vaginitis_status) & 
  metadata_all_2$hepatitisB_status == 'No' & !is.na(metadata_all_2$hepatitisB_status) & 
  metadata_all_2$hepatitisC_status == 'No' & !is.na(metadata_all_2$hepatitisC_status) & 
  metadata_all_2$hiv_status == 'No' & !is.na(metadata_all_2$hiv_status) & 
  metadata_all_2$cancer_status == 'No' & !is.na(metadata_all_2$cancer_status) & 
  metadata_all_2$douche == 0 & !is.na(metadata_all_2$douche)


keep = metadata_all_2$quinolones_status == 'No' & !is.na(metadata_all_2$quinolones_status) & 
  metadata_all_2$sulfa_status == 'No' & !is.na(metadata_all_2$sulfa_status) & 
  metadata_all_2$penicillins_status == 'No' & !is.na(metadata_all_2$penicillins_status) & 
  metadata_all_2$other_antibiotics_status == 'No' & !is.na(metadata_all_2$other_antibiotics_status) & 
  metadata_all_2$cephalosporins_status == 'No' & !is.na(metadata_all_2$cephalosporins_status) & 
  metadata_all_2$macrolides_status == 'No' & !is.na(metadata_all_2$macrolides_status) & 
  metadata_all_2$nitrofurantoin_status == 'No' & !is.na(metadata_all_2$nitrofurantoin_status) & 
  metadata_all_2$insulin_status == 'No' & !is.na(metadata_all_2$insulin_status) & 
  metadata_all_2$oral_steroids_status == 'No' & !is.na(metadata_all_2$oral_steroids_status) & 
  metadata_all_2$tetracyclines_status == 'No' & !is.na(metadata_all_2$tetracyclines_status) & 
  metadata_all_2$metrodinazole_status == 'No' & !is.na(metadata_all_2$metrodinazole_status) & 
  metadata_all_2$hiv_status == 'No' & !is.na(metadata_all_2$hiv_status) & 
  metadata_all_2$cancer_status == 'No' & !is.na(metadata_all_2$cancer_status) 

sum(keep)

reads_table_all_abundance_3 = reads_table_all_abundance_2[,keep]
metadata_all_3 = metadata_all_2[keep,]
mytypes_3 = mytypes[keep]

reads_table = as.data.frame(t(reads_table_all_abundance_3))
reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
mytypes_2 <- apply(reads_table,1,which.max)   # find the most abundant taxa
maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes_2), ncol=2)]  # the abundance of the most abundant taxa

sampleID_list = as.data.frame(matrix(data = NA, ncol = 14, nrow = 8))
colnames(sampleID_list) = c(top_10,'G.vag_bv')

for (a in 1: 13) {
  keep = mytypes_3 == top_10[a]
  if (sum(keep) < 8) {
    print(paste0(top_10[a],'    ',sum(keep)))
    break
  }
  
  maxprop_2 = maxprop[keep]
  sampleID_list_2 = row.names(metadata_all_3)[keep]
  
  if (top_10[a] != "No_Type") {
    n = order(maxprop_2, decreasing = T)
  } else {
    n = order(maxprop_2, decreasing = F)
  }
  maxprop_2[n[1:8]]
  sampleID_list_2 = sampleID_list_2[n[1:8]]
  sampleID_list[,a] = sampleID_list_2
}
sampleID_list$Prevotella_cluster2[c(3:8)] = NA
sampleID_list$Gardnerella_vaginalis = NA

# G. vag
keep = metadata_all_2$quinolones_status == 'No' & !is.na(metadata_all_2$quinolones_status) & 
  metadata_all_2$sulfa_status == 'No' & !is.na(metadata_all_2$sulfa_status) & 
  metadata_all_2$penicillins_status == 'No' & !is.na(metadata_all_2$penicillins_status) & 
  metadata_all_2$other_antibiotics_status == 'No' & !is.na(metadata_all_2$other_antibiotics_status) & 
  metadata_all_2$cephalosporins_status == 'No' & !is.na(metadata_all_2$cephalosporins_status) & 
  metadata_all_2$macrolides_status == 'No' & !is.na(metadata_all_2$macrolides_status) & 
  metadata_all_2$nitrofurantoin_status == 'No' & !is.na(metadata_all_2$nitrofurantoin_status) & 
  metadata_all_2$insulin_status == 'No' & !is.na(metadata_all_2$insulin_status) & 
  metadata_all_2$oral_steroids_status == 'No' & !is.na(metadata_all_2$oral_steroids_status) & 
  metadata_all_2$tetracyclines_status == 'No' & !is.na(metadata_all_2$tetracyclines_status) & 
  metadata_all_2$metrodinazole_status == 'No' & !is.na(metadata_all_2$metrodinazole_status) & 
  metadata_all_2$hiv_status == 'No' & !is.na(metadata_all_2$hiv_status) & 
  metadata_all_2$cancer_status == 'No' & !is.na(metadata_all_2$cancer_status) & 
  metadata_all_2$dx_bv == 'No' & !is.na(metadata_all_2$dx_bv)& 
  metadata_all_2$bv == 'No' & !is.na(metadata_all_2$bv)& 
  metadata_all_2$bacterial_vaginosis_status == 'No' & !is.na(metadata_all_2$bacterial_vaginosis_status)

sum(keep)

reads_table_all_abundance_3 = reads_table_all_abundance_2[,keep]
metadata_all_3 = metadata_all_2[keep,]
mytypes_3 = mytypes[keep]

reads_table = as.data.frame(t(reads_table_all_abundance_3))
reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
mytypes_2 <- apply(reads_table,1,which.max)   # find the most abundant taxa
maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes_2), ncol=2)]  # the abundance of the most abundant taxa

a=5
  keep = mytypes_3 == top_10[a]
  if (sum(keep) < 8) {
    print(paste0(top_10[a],'    ',sum(keep)))
    break
  }
  
  maxprop_2 = maxprop[keep]
  sampleID_list_2 = row.names(metadata_all_3)[keep]
  
  if (top_10[a] != "No_Type") {
    n = order(maxprop_2, decreasing = T)
  } else {
    n = order(maxprop_2, decreasing = F)
  }
  maxprop_2[n[1:8]]
  sampleID_list_2 = sampleID_list_2[n[1:8]]
  sampleID_list[,a] = sampleID_list_2

# G. vag bv
  keep = metadata_all_2$quinolones_status == 'No' & !is.na(metadata_all_2$quinolones_status) & 
    metadata_all_2$sulfa_status == 'No' & !is.na(metadata_all_2$sulfa_status) & 
    metadata_all_2$penicillins_status == 'No' & !is.na(metadata_all_2$penicillins_status) & 
    metadata_all_2$other_antibiotics_status == 'No' & !is.na(metadata_all_2$other_antibiotics_status) & 
    metadata_all_2$cephalosporins_status == 'No' & !is.na(metadata_all_2$cephalosporins_status) & 
    metadata_all_2$macrolides_status == 'No' & !is.na(metadata_all_2$macrolides_status) & 
    metadata_all_2$nitrofurantoin_status == 'No' & !is.na(metadata_all_2$nitrofurantoin_status) & 
    metadata_all_2$insulin_status == 'No' & !is.na(metadata_all_2$insulin_status) & 
    metadata_all_2$oral_steroids_status == 'No' & !is.na(metadata_all_2$oral_steroids_status) & 
    metadata_all_2$tetracyclines_status == 'No' & !is.na(metadata_all_2$tetracyclines_status) & 
    metadata_all_2$metrodinazole_status == 'No' & !is.na(metadata_all_2$metrodinazole_status) & 
    metadata_all_2$hiv_status == 'No' & !is.na(metadata_all_2$hiv_status) & 
    metadata_all_2$cancer_status == 'No' & !is.na(metadata_all_2$cancer_status) & 
    metadata_all_2$dx_bv == 'Yes' & !is.na(metadata_all_2$dx_bv)& 
    metadata_all_2$bv == 'Yes' & !is.na(metadata_all_2$bv)& 
    metadata_all_2$bacterial_vaginosis_status == 'Yes' & !is.na(metadata_all_2$bacterial_vaginosis_status)
  
  sum(keep)
  
  reads_table_all_abundance_3 = reads_table_all_abundance_2[,keep]
  metadata_all_3 = metadata_all_2[keep,]
  mytypes_3 = mytypes[keep]
  
  reads_table = as.data.frame(t(reads_table_all_abundance_3))
  reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
  mytypes_2 <- apply(reads_table,1,which.max)   # find the most abundant taxa
  maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes_2), ncol=2)]  # the abundance of the most abundant taxa
  
  a=5
  keep = mytypes_3 == top_10[a]
  if (sum(keep) < 8) {
    print(paste0(top_10[a],'    ',sum(keep)))
    break
  }
  
  maxprop_2 = maxprop[keep]
  sampleID_list_2 = row.names(metadata_all_3)[keep]
  
  if (top_10[a] != "No_Type") {
    n = order(maxprop_2, decreasing = T)
  } else {
    n = order(maxprop_2, decreasing = F)
  }
  maxprop_2[n[1:8]]
  sampleID_list_2 = sampleID_list_2[n[1:8]]
  sampleID_list[,14] = sampleID_list_2



setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu')
write.csv(sampleID_list,'VAMP_10Sample_for_top10_vagitype.csv', quote = F)

setwd('/Users/binzhu/Desktop/pH/results')
sampleID_list_2 = gather(sampleID_list)
sampleID_list_2 = sampleID_list_2$value

setwd('/Users/binzhu/secure/fenn/vcu_gpfs2/home/bzhu/bccl/bzhu')
write.table(sampleID_list_2,'VAMP_10Sample_for_top10_vagitype.txt', quote = F, row.names = F, col.names = F)

########## others ##########
n = row.names(reads_table_all) %in% taxa_order$Taxa
sum(rowSums(reads_table_all)[n]) / sum(rowSums(reads_table_all))
sum(BVT != 'Others') / length(BVT)
