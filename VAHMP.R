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
#reads_table = reads_table$reads_table
row.names(reads_table)[row.names(reads_table) == 'Lactobacillus_crispatus_cluster'] = 'Lactobacillus_crispatus'
row.names(reads_table)[row.names(reads_table) == 'Lactobacillus_gasseri_cluster'] = 'Lactobacillus_gasseri'
row.names(reads_table)[row.names(reads_table) == 'Lachnospiraceae_BVAB1'] = 'Ca._L._vaginae'
row.names(reads_table)[row.names(reads_table) == 'BVAB'] = 'BVAB_TM7'
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
#reads_table_all = reads_table_all[-which(row.names(reads_table_all) == 'Other_Prevotellaceae'),]

#write.csv(row.names(reads_table_all),'taxa_list.csv')
#write.csv(reads_table_all,'reads_table_all.csv')
#write.csv(metadata_all,'metadata_all.csv')









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
   
   {
     Median_relative_abundance = as.data.frame(mytypes)
     Median_relative_abundance$relative_abundance = NA
     
     for (a in 1:nrow(Median_relative_abundance)) {
       n = which(row.names(reads_table_all_abundance) == as.character(Median_relative_abundance$mytypes[a]))
       
       if (length(n) == 0) {
         next
       }
       Median_relative_abundance$relative_abundance[a] = reads_table_all_abundance[n,a]
     }
     Median_relative_abundance = Median_relative_abundance[Median_relative_abundance$mytypes != 'Others' & Median_relative_abundance$mytypes != 'No_Type',]
     
     Median_relative_abundance_taxa_list = unique(Median_relative_abundance$mytypes)
     Median_relative_abundance_taxa_list = as.data.frame(Median_relative_abundance_taxa_list)
     Median_relative_abundance_taxa_list$Median = NA
     
     for (a in 1: nrow(Median_relative_abundance_taxa_list)) {
       n = which(Median_relative_abundance$mytypes == Median_relative_abundance_taxa_list$Median_relative_abundance_taxa_list[a])
       median(Median_relative_abundance$relative_abundance[n])
       Median_relative_abundance_taxa_list$Median[a] = median(Median_relative_abundance$relative_abundance[n])
     }
     write.csv(Median_relative_abundance_taxa_list,'Median_relative_abundance_taxa_list.csv')
   }
   
   data = cbind(colnames(reads_table_all), mytypes)
   write.csv(data,'Vagitype.csv')
   top_vagitype = as.data.frame(table(mytypes))
   top_vagitype = top_vagitype[order(top_vagitype$Freq, decreasing = T),]
   
   write.csv(top_vagitype,'top_vagitype.csv')
   keep = mytypes %in% top_vagitype$mytypes[1:12]
   
   mytypes[!keep] = 'Others'
   metadata_all$Vagitype = mytypes
   
   top_vagitype$Median_relative_abundance = NA
   
   Median_relative_abundance = as.data.frame(mytypes)
   Median_relative_abundance$relative_abundance = NA
   
   for (a in 1:nrow(Median_relative_abundance)) {
     n = which(row.names(reads_table_all_abundance) == as.character(Median_relative_abundance$mytypes[a]))
     
     if (length(n) == 0) {
       next
     }
     Median_relative_abundance$relative_abundance[a] = reads_table_all_abundance[n,a]
   }
   Median_relative_abundance = Median_relative_abundance[Median_relative_abundance$mytypes != 'Others' & Median_relative_abundance$mytypes != 'No_Type',]
   
   Median_relative_abundance$mytypes = str_replace_all(Median_relative_abundance$mytypes, "_" ,' ')
   Median_relative_abundance$mytypes = factor(Median_relative_abundance$mytypes, levels = 
                                                   c("Lactobacillus crispatus","Lactobacillus gasseri","Lactobacillus iners","Lactobacillus jensenii",
                                                     "Atopobium vaginae","Gardnerella vaginalis","Ca. L. vaginae", "Streptococcus", "Mycoplasma","Sneathia amnii", "Prevotella bivia"))

   ggplot(Median_relative_abundance, aes(x=mytypes, y=relative_abundance, fill = mytypes)) + 
     geom_boxplot(color="black", outlier.shape=NA, width=0.5) +
     scale_fill_manual(values=c("#FFDC00","#c5b0d5", 
                                "#aec7e8","#49C800","#c49c94","#d62728", "#ff7f0e", 
                                "#7f7f7f","#e377c2", "#FF00D1", "#17becf"))+
     labs(x ='Vagitype', y = "Relative abundance")+ theme_bw()+
     theme(axis.title = element_text(size = 7), 
           axis.text = element_text(size = 7), 
           legend.text = element_text(size = 7), 
           legend.title = element_text(size = 7),
           axis.text.x = element_text(angle = 65, vjust = 1, hjust=1))
   ggsave('Median_relative_abundance_when_abundant.pdf',width = 4,height = 2.5)
   
   taxa_list = unique(Median_relative_abundance$mytypes)
   output = matrix(data = NA, nrow = length(taxa_list), ncol = length(taxa_list))
   colnames(output) = taxa_list; row.names(output) = taxa_list
   {
     for (a in 1: (length(taxa_list)-1)) {
       for (b in (a+1) : length(taxa_list)) {
         data = Median_relative_abundance[Median_relative_abundance$mytypes == taxa_list[a] | Median_relative_abundance$mytypes == taxa_list[b],]
         x = wilcox.test(relative_abundance~mytypes, data = data)
         output[a,b] = x$p.value
       }
     }
   }
   write.csv(output,'Median_relative_abundance_when_abundant.csv')
   
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
    group_level = as.character(group_level)
    
    within_dis_sig <- matrix(data = NA, ncol=n, nrow = n)
    colnames(within_dis_sig) = group_level
    row.names(within_dis_sig) = group_level
    for (a in 1:(n-1)) {
      for (b in (a+1): n) {
        
        keep = metadata == group_level[a] | metadata == group_level[b]
        metadata_2 = metadata[keep]
        reads_table_2 = reads_table[keep,]
        data = mrpp(reads_table_2, as.matrix(metadata_2), permutations = 999, distance = "bray",
                    weight.type = 1, strata = NULL, parallel = getOption("mc.cores"))
        within_dis_sig[a,b] <- data$Pvalue
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

###### correlation between taxa and pH #####
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

  p <- list(); output = data.frame(Taxa = row.names(reads_table), R = NA, P = NA, adj_p = NA)
  
  for (a in seq(nrow(reads_table))) {
    data = data.frame(pH = metadata$sample_ph, relative_abundance = as.numeric(as.character(reads_table[a,])))
    
    # linear
    linearMod = lm(relative_abundance~pH, data)
    
    pvalue <- as.numeric(summary(linearMod)$coefficients[,4][2])  # get p-value
    output$P[a] = pvalue
    line_list[a,2]= pvalue
    r_squared <- summary(linearMod)$r.squared  # get R-value
    output$R[a] = r_squared
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
  data <- adjust.p(output$P, pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
  data <- data$adjp
  data <- data$adjusted.p
  output$adj_p = data
  
  for (a in seq(nrow(reads_table))) {
    data = data.frame(pH = metadata$sample_ph, relative_abundance = as.numeric(as.character(reads_table[a,])))
    
    # linear
    linearMod = lm(relative_abundance~pH, data)
    
    pvalue <- output$adj_p[a] 
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
      labs(title = paste0(str_replace_all(row.names(reads_table)[a],"_"," "),'\n',line_list[a,4],'\n',"Adj-P-value = ",formatC(pvalue, format = "e", digits = 2),
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
  paletteLength <- 50
  pvalue = 0.05;cor_parameter= 0;style = 1; pheatmap_fontsize = 5
  treeheight = 50; alpha = 0.05; mc.cores =60
  bar_max = 0.6
  bar_min = -0.6
  
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
                     treeheight_row = 25, treeheight_col = 25, method = "complete",
                     clustering_distance_rows = "euclidean",
                     clustering_distance_cols = "euclidean",show_colnames = T, 
                     annotation_colors = mycolors,annotation_row = taxa_cluster,annotation_col = taxa_cluster)
  
  pdf("network_all.pdf", width=6, height=5)
  print(cluster)
  dev.off()
  
  x = dist(p.yes.rr, method = "euclidean", diag = FALSE, upper = T)
  x= as.matrix(x)
  write.csv(x, 'distance_in_network.csv')
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
  
  paletteLength <- 50
  pvalue = 0.05;cor_parameter= 0;style = 1; pheatmap_fontsize = 5
  treeheight = 50; alpha = 0.05; mc.cores =60
  bar_max = max(p.yes.rr,na.rm = T)
  bar_min = min(p.yes.rr,na.rm = T)
  
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

{
  data = read.csv('model_pH_network_MTS.csv',row.names = 1)
  data = data[!is.na(data$Sensitivity_to_pH),]
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
  pdf("Taxa_classification_2.pdf", width=12, height=4.5)
  old.par <- par(mfrow=c(1, 2))
  plot(hcl, labels=gsub("[A-z.]", "  ", hcl$labels))
  Tctext(hcl, srt=90, add=0.04, adj=c(1, 0.5))
  plot(hcl)
  par(old.par)
  dev.off()
  
  pdf("Taxa_classification_1.pdf", width=12, height=8)
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
  
  mycolors = list(Cluster = c("I"="#FFDC00","II"="#49C800","III"="#aec7e8","IV"="#17becf","V"="#7f7f7f",
                              "VI"="#ff7f0e","VII"="#d62728","Others"="#9C9146"))
  
  cluster = pheatmap(data_2, cluster_rows=F, show_rownames=TRUE,
                     cluster_cols=F, annotation_row = taxa_order,annotation_colors =mycolors)
  
  pdf("Taxa_classification_3.pdf", width=3.3, height=5)
  print(cluster)
  dev.off()
  
}


########## application of B-vagitype ############
setwd('/Users/binzhu/Desktop/pH/results')
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
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 9] = 'VI'
  
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
                          "VI"="#e377c2"),
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

# prepare metadata
{
  
  colnames(metadata_all)[colnames(metadata_all) == 'mom_delivery_method'] = 'mom_delivery_method (C-section yes / vaginal no)'
  colnames(metadata_all)[colnames(metadata_all) == 'Delivery'] = 'Preterm delivery'
  colnames(metadata_all)[colnames(metadata_all) == 'baby_icu'] = 'baby_nicu'
  colnames(metadata_all)[colnames(metadata_all) == 'baby_problems_none'] = 'baby_problems'
  
  metadata_all[metadata_all == 'PTB' | metadata_all == 'Yes_I_plan_to_breast_feed' | 
                 metadata_all == 'Yes_my_baby_was_put_in_the_ICU' |metadata_all == 'Yes_I_have_been_diagnosed_with_hyperemesis' |
                 metadata_all == 'C_Section' |metadata_all == 'Cesarean_section' |  
                 metadata_all == 'Female'| metadata_all == 'Yes_98_100'| metadata_all == 'Yes_101_103'| metadata_all == 'Yes_membranes_ruptured_1_18_hours_before_labor'| metadata_all == 'Yes_membranes_ruptured_before_37_weeks_of_pregnancy'] = 'Yes'
  
  metadata_all[metadata_all == 'TB' | metadata_all == 'No_I_have_not_been_exposed_to_any_chemicals' | metadata_all == 'No_I_have_not_been_to_the_hospital' | metadata_all == 'No_I_have_not_been_diagnosed_with_hyperemesis_and_I_do_not_have_symptoms' |
                 metadata_all == 'No_my_baby_was_not_put_in_the_ICU' | metadata_all == 'No_I_do_not_plan_to_breast_feed' | metadata_all == 'No_I_am_not_taking_any_new_medications' | metadata_all == 'No_I_have_not_been_diagnosed_with_hyperemesis_but_I_have_symptoms' |
                 metadata_all == 'Vaginal' | metadata_all == 'Vaginal_delivery' |  metadata_all == 'No_I_have_not_had_a_fever_since_my_last_visit' |
                 metadata_all == 'Male' | metadata_all == 'none'| metadata_all == 'No_there_has_been_no_premature_rupture_of_membranes'] = 'No'
  
  metadata_all[metadata_all == 'Not_Sure' | 
                 metadata_all == 'Not_sure' | 
                 metadata_all == 'Uncertain' |
                 metadata_all == 'Do_not_know'|
                 metadata_all == 'Nelson 6th Floor'|
                 metadata_all == 'Baby_NICU' | metadata_all == 'NaN' | metadata_all == 'Don_t_know'] = NA
  
  # write.csv(metadata_all,'metadata_all_ori.csv') ### run when output contains sample ID or participant ID
  
  # smoker_second_hand
  metadata_all[metadata_all == 'Never'] = 0
  metadata_all[metadata_all == 'Rarely'] = 1
  metadata_all[metadata_all == 'Almost_every_day'] = 2
  metadata_all[metadata_all == 'Every_day'] = 3
  
  # various worries
  metadata_all[metadata_all == 'Much_Less'] = 1
  metadata_all[metadata_all == 'A_Bit_Less'] = 2
  metadata_all[metadata_all == 'About_Average'] = 3
  metadata_all[metadata_all == 'A_Bit_More'] = 4
  metadata_all[metadata_all == 'Much_More'] = 5
  
  metadata_all[metadata_all == 'Never'] = 0
  metadata_all[metadata_all == 'Almost_Never'] = 1
  metadata_all[metadata_all == 'Sometimes'] = 2
  metadata_all[metadata_all == 'Fairly_Often'] = 3
  metadata_all[metadata_all == 'Very_Often'] = 4
  
  
  # yogurt, milk, cheese, ice_cream
  metadata_all[metadata_all == '1_2_servings_per_day'] = 4
  metadata_all[metadata_all == '3_4_servings_per_day'] = 5
  metadata_all[metadata_all == '5_servings_per_day'] = 6
  metadata_all[metadata_all == 'Less_than_1_serving_per_week'] = 1
  metadata_all[metadata_all == '1_2_servings_per_week'] = 2
  metadata_all[metadata_all == '3_4_servings_per_week'] = 3
  
  # weight_change
  metadata_all[metadata_all == 'My_weight_has_NOT_changed_since_my_last_visit'] = 0
  metadata_all[metadata_all == 'I_have_LOST_weight_since_my_last_visit'] = -1
  metadata_all[metadata_all == 'I_have_GAINED_weight_since_my_last_visit'] = 1
  
  # weight_lost, weight_gained
  metadata_all[metadata_all == '0_5_lbs'] = 1
  metadata_all[metadata_all == '10_15_lbs'] = 3
  metadata_all[metadata_all == '5_10_lbs'] = 2
  metadata_all[metadata_all == 'More_than_15_lbs'] = 4
  
  # physical_activity_vigorous, physical_activity_moderate, physical_activity_light
  metadata_all[metadata_all == '0_times'] = 0
  metadata_all[metadata_all == '1_2_times'] = 1
  metadata_all[metadata_all == '3_4_times'] = 2
  metadata_all[metadata_all == '5_6_times'] = 3
  metadata_all[metadata_all == '7_times'] = 4
  
  # uti_lifetime_number, 
  metadata_all[metadata_all == '2_to_4'] = 3
  metadata_all[metadata_all == '2_to_4'] = 3
  
  # douche_frequency, vaginal_sex_frequency, vaginal_penetration_frequency, received_oral_frequency, gave_oral_frequency, anal_sex_frequency, vaginal_lubrication_frequency
  metadata_all[metadata_all == '1_3_times_per_month'] = 2
  metadata_all[metadata_all == '2_6_times_per_week'] = 4
  metadata_all[metadata_all == 'Less_than_once_a_month'] = 1
  metadata_all[metadata_all == 'Never'] = 0
  metadata_all[metadata_all == 'Once_a_day'] = 5
  metadata_all[metadata_all == 'Once_a_week'] = 3
  
  # income
  metadata_all[metadata_all == '15_000_19_999'] = 2
  metadata_all[metadata_all == '20_000_39_999'] = 3
  metadata_all[metadata_all == '40_000_59_999'] = 4
  metadata_all[metadata_all == '60_000_79_999'] = 5
  metadata_all[metadata_all == '80_000_or_more'] = 6
  metadata_all[metadata_all == 'Under_15_000'] = 1
  
  # education
  metadata_all[metadata_all == '2_year_College_Degree'] = 3
  metadata_all[metadata_all == '4_year_College_Degree'] = 4
  metadata_all[metadata_all == 'Doctoral_or_Professional_Degree'] = 6
  metadata_all[metadata_all == 'High_School_GED'] = 2
  metadata_all[metadata_all == 'Less_than_High_School'] = 1
  metadata_all[metadata_all == 'Masters_Degree'] = 5
  metadata_all[metadata_all == 'Some_College'] = NA
  
  # prenatal_care_start
  metadata_all[metadata_all == '0_4_weeks_after_conception'] = 3
  metadata_all[metadata_all == '4_weeks_to_the_end_of_the_1st_trimester'] = 2
  metadata_all[metadata_all == '4_weeks_to_the_end_of_the_first_trimester'] = 2
  metadata_all[metadata_all == 'Before_conception'] = 4
  metadata_all[metadata_all == 'During_the_2nd_trimester'] = 1
  metadata_all[metadata_all == 'During_the_second_trimester'] = 1
  
  metadata_all[metadata_all == 'Yes'] =1
  metadata_all[metadata_all == 'No'] =0
  
  metadata = metadata_all
  metadata = apply(metadata, 2, as.numeric)
  keep = sapply(1:ncol(metadata), function(j) (length(unique(metadata[,j])) >= 2 ))
  metadata = metadata[,keep]
  
  metadata = as.data.frame(metadata)
  
  for (a in 1: ncol(metadata)) {
    data = metadata[,a]
    keep = unique(data); keep = keep[!is.na(keep)]
    if (length(keep) >2 & length(keep) <=6) {
      keep = data > mean(data[!is.na(data)]) & !is.na(data); metadata[keep,a] = 1
      keep = data <= mean(data[!is.na(data)]) & !is.na(data); metadata[keep,a] = 0
    } else {next}
    
  }
  
}

# test classification
{
  # modulation of variables
  {
    # find correlated variables
    metadata_list = read.csv('metadata_list.csv',header = F)
    metadata_list = metadata_list$V1
    metadata = metadata[,colnames(metadata) %in% metadata_list]
    data =  t(metadata)
    output = newwork_rcorr(data, normalization_method = NA, type = 'spearman', pvalue = 0.05, 
                           cor_parameter= 0, style = 1, bar_max = 1, bar_min = -1, 
                           pheatmap_fontsize = 5, treeheight = 20, alpha = 0.05)
    p.yes.rr = abs(output$cor_matrix)
    
    bar_max = max(p.yes.rr,na.rm = T)
    bar_min = min(p.yes.rr,na.rm = T)
    paletteLength <- 50
    myColor <- colorRampPalette(c("white", "red"))(paletteLength)
    p = pheatmap(p.yes.rr, color=myColor, fontsize = 5, 
                 treeheight_row = 50, treeheight_col = 50)
    
    pdf("correlated_variables.pdf", width=12, height=12)
    print(p)
    dev.off()
    
    correlated_variables_network = output$gephi_input
    correlated_variables_network$Weight = abs(correlated_variables_network$Weight)
    correlated_variables_network$Weight = as.numeric(as.character(correlated_variables_network$Weight))
    write.csv(correlated_variables_network,'correlated_variables_network.csv', row.names = F)
    
    data_3 = as.data.frame(output$cor_matrix)
    data_3 = gather(data_3)
    data_4 = as.data.frame(output$adj_p)
    data_4 = gather(data_4)
    data_3$adj_p = data_4$value
    data_3$variable_2 = rep(metadata_list,ncol(output$cor_matrix))
    colnames(data_3) = c('Variable 1', 'R-value','Adj-P-value','Variable 2')
    data_3 = data_3[!is.na(data_3$`Adj-P-value`),]
    write.csv(data_3,'correlated_variables_network_2.csv',row.names = F)
    
    p_cluster = read.csv('gephi_metadata_network.csv')
    #    p_cluster = sort(cutree(p$tree_col, k=30))
    #    p_cluster = as.data.frame(p_cluster)
    
    p_cluster = p_cluster[order(p_cluster$betweenesscentrality,decreasing = T),]


  }
  # test classification
  metadata_test = data.frame(CST = CST, Vagitype= mytypes, BVT = BVT) 
  for (n in 1:3) {
    output_table = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = length(unique(metadata_test[,n]))))
    
    colnames(output_table) = colnames(metadata)
    row.names(output_table) = sort(unique(metadata_test[,n]))
    
    output_table_2 = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = 2))
    colnames(output_table_2) = colnames(metadata)
    row.names(output_table_2) = c('Kruskal–Wallis test',"Pearson's Chi-squared test")
    
    output_table_3 = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = 1))
    colnames(output_table_3) = colnames(metadata)
    
    for (y in 1: ncol(output_table)) {
      data = metadata[,y]
      data = data[!is.na(data)]
      if (length(unique(data)) > 6) {
        data_2 = median(data)
        data_2 = round(data_2, 1)
        data_3 = quantile(data)
        data_3 = round(data_3, 1)
        output_table_3[1,y] = paste0(data_2, " (",data_3[2], ", ",data_3[4], "); ",length(data))
      } else if (length(unique(data)) == 2) {
        data_2 = sum(data == 1) / length(data)
        data_2 = round(data_2, 3)
        output_table_3[1,y] = paste0(data_2*100,'%; ',length(data))
      } else {next}
      
      for (x in 1:nrow(output_table)) {
        data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
        data = data[!is.na(data)]
        
        if (length(unique(data)) > 6) {
          data_2 = median(data)
          data_2 = round(data_2, 1)
          data_3 = quantile(data)
          data_3 = round(data_3, 1)
          output_table[x,y] = paste0(data_2, " (",data_3[2], ", ",data_3[4], "); ",length(data))
        } else if (length(unique(data)) == 2) {
          data_2 = sum(data == 1) / length(data)
          data_2 = round(data_2, 3)
          output_table[x,y] = paste0(data_2*100,'%; ',length(data))
        } else {next}
      }
      
      data_4 = as.data.frame(cbind(metadata_test[,n], metadata[,y]))
      keep = !is.na(data_4[,2])
      data_4 = data_4[keep,]
      
      if (length(unique(data_4$V2)) == 2 | length(unique(data_4$V2)) >6) {
        output_table_2[1,y] = kruskal.test(V2~V1, data_4)$p.value
        output_table_2[2,y] = chisq.test(data_4$V2,  data_4$V1, correct=FALSE)$p.value
      }
      
    }
    write.csv(output_table_3,'Sample_size.csv',row.names = F)
    
    output_table = rbind(output_table,output_table_2)
    
    if (n ==1) {output_table_CST = output_table} else if (n ==2) {output_table_vagitype = output_table}  else {output_table_BVT = output_table} 
  }
  {
    keep =  colSums(is.na(output_table_CST)) < nrow(output_table_CST)/3 & 
      colSums(is.na(output_table_vagitype)) < nrow(output_table_vagitype)/3 & 
      colSums(is.na(output_table_BVT)) < nrow(output_table_BVT)/3
    
    sum(keep)
    
    output_table_CST = output_table_CST[ , keep]
    output_table_vagitype = output_table_vagitype[ ,keep]
    output_table_BVT = output_table_BVT[ ,keep]
    
    compare_list = data.frame(Metadata = colnames(output_table_CST), CST_pvalue = as.numeric(as.character(output_table_CST[11,])),
                              Vagitype_pvalue = as.numeric(as.character(output_table_vagitype[15,])), 
                              BVT_pvalue = as.numeric(as.character(output_table_BVT[10,])), Best_pvalue = NA)
    
    for (a in 2:4) {
      data <- adjust.p(compare_list[,a], pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data <- data$adjp
      data <- data$adjusted.p
      compare_list[,a] = data
    }
    output_table_CST[11,] = compare_list[,2];output_table_vagitype[15,] = compare_list[,3];output_table_BVT[10,] = compare_list[,4]
    
    compare_list_3 = data.frame(Metadata = colnames(output_table_CST), CST_pvalue = as.numeric(as.character(output_table_CST[10,])),
                                Vagitype_pvalue = as.numeric(as.character(output_table_vagitype[14,])), 
                                BVT_pvalue = as.numeric(as.character(output_table_BVT[9,])), Best_pvalue = NA)
    
    for (a in 2:4) {
      data <- adjust.p(compare_list_3[,a], pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data <- data$adjp
      data <- data$adjusted.p
      compare_list_3[,a] = data
    }
    output_table_CST[10,] = compare_list_3[,2];output_table_vagitype[14,] = compare_list_3[,3];output_table_BVT[9,] = compare_list_3[,4]
    
    
    compare_list_2 = sapply(1:nrow(compare_list), function(j) (compare_list[j,5] = which.min(compare_list[j,c(2:4)])))
    compare_list_2[compare_list_2 == 1] = 'CST';compare_list_2[compare_list_2 == 2] = 'Vagitype';compare_list_2[compare_list_2 == 3] = 'BVT';
    compare_list[,5] = compare_list_2
    
    keep = compare_list[,2] <= 0.05 | compare_list[,3] <= 0.05 | compare_list[,4] <= 0.05; sum(keep)
    compare_list[!keep,5] = ''
    
    metadata_list = compare_list$Metadata[!is.na(compare_list$Best_pvalue)]
    #  metadata_list = read.csv('metadata_list.csv', header = F)
    #  metadata_list = metadata_list$V1
    output_table_CST = output_table_CST[,colnames(output_table_CST) %in% metadata_list]
    output_table_vagitype = output_table_vagitype[,colnames(output_table_vagitype) %in% metadata_list]
    output_table_BVT = output_table_BVT[,colnames(output_table_BVT) %in% metadata_list]
    #  compare_list = compare_list[compare_list$Metadata %in% metadata_list,]
    
    write.csv(output_table_CST, paste0('Distribution_CST.csv'))
    write.csv(output_table_vagitype, paste0('Distribution_Vagitype.csv'))
    write.csv(output_table_BVT, paste0('Distribution_BVT.csv'))
    output_table = rbind(output_table_CST,output_table_vagitype,output_table_BVT)
    write.csv(output_table, paste0('Distribution_all.csv'))
    write.csv(compare_list, paste0('compare_list.csv'))
    
    data = compare_list[,c(2:5)]
    data = gather(data)
    data$key = str_remove_all(data$key,'_pvalue')
    colnames(data) = c('Method','-Log10(Adj-P-value)')
    data$`-Log10(Adj-P-value)` = -log10(as.numeric(as.character(data$`-Log10(Adj-P-value)`)))
    data = data[!is.na(data$`-Log10(Adj-P-value)`),]
    
    data$Method = factor(data$Method, levels = c('CST', 'Vagitype', 'BVT'))
    ggplot(data, aes(x = Method, y = `-Log10(Adj-P-value)`, color = Method)) +
      geom_boxplot(outlier.shape = NA)+
      geom_jitter(size = 0.1)+ theme_bw()+
      theme(axis.title = element_text(size = 7), 
            axis.text = element_text(size = 7), 
            legend.text = element_text(size = 7), 
            legend.title = element_text(size = 7))
    ggsave('Method.compare.pdf',height = 3, width = 3)
    
    # calculate significance
    {
      data$Method = as.character(data$Method)
      factor_levels = unique(data$Method)
      n = length(factor_levels)
      
      Method_sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
      colnames(Method_sig) = factor_levels
      row.names(Method_sig) = factor_levels
      
      for (a in 1:(n-1)) {
        for (b in (a+1) : n) {
          factor_level1 <- data$`-Log10(Adj-P-value)`[data$Method == factor_levels[a]]
          factor_level2 <- data$`-Log10(Adj-P-value)`[data$Method == factor_levels[b]]
          
          Method_sig[a,b] <- wilcox.test(factor_level1, factor_level2, paired = T)$p.value
          
        }
      }
      write.csv(Method_sig,'Method.compare.csv')
    }
    
    # krus test
    {
      data = compare_list_3[,c(2:5)]
      data = gather(data)
      data$key = str_remove_all(data$key,'_pvalue')
      colnames(data) = c('Method','-Log10(Adj-P-value)')
      data$`-Log10(Adj-P-value)` = -log10(as.numeric(as.character(data$`-Log10(Adj-P-value)`)))
      data = data[!is.na(data$`-Log10(Adj-P-value)`),]
      
      data$Method = factor(data$Method, levels = c('CST', 'Vagitype', 'BVT'))
      ggplot(data, aes(x = Method, y = `-Log10(Adj-P-value)`, color = Method)) +
        geom_boxplot(outlier.shape = NA)+
        geom_jitter(size = 0.1)+ theme_bw()+
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))
      ggsave('Method.compare_krus.pdf',height = 3, width = 3)
      
      # calculate significance
      {
        data$Method = as.character(data$Method)
        factor_levels = unique(data$Method)
        n = length(factor_levels)
        
        Method_sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
        colnames(Method_sig) = factor_levels
        row.names(Method_sig) = factor_levels
        
        for (a in 1:(n-1)) {
          for (b in (a+1) : n) {
            factor_level1 <- data$`-Log10(Adj-P-value)`[data$Method == factor_levels[a]]
            factor_level2 <- data$`-Log10(Adj-P-value)`[data$Method == factor_levels[b]]
            
            Method_sig[a,b] <- wilcox.test(factor_level1, factor_level2, paired = T)$p.value
            
          }
        }
        write.csv(Method_sig,'Method.compare_krus.csv')
      }
    }
    
  }
  p_cluster$modularity_class = as.factor(p_cluster$modularity_class)
  p_cluster$BVT_adj_p_pearson = compare_list$BVT_pvalue[match(p_cluster$Label, compare_list$Metadata)]
  p_cluster$BVT_adj_p_pearson = -log10(p_cluster$BVT_adj_p_pearson)
  ggplot(p_cluster, aes(x = modularity_class, y =BVT_adj_p_pearson, color = modularity_class)) +
    geom_boxplot(outlier.shape = NA)+geom_hline(yintercept=-log10(0.05),linetype=2)+
    geom_jitter(size = 0.1)+ theme_bw()+ 
    scale_color_manual(values=c("#2cb384", "#34b5ff", "#38d011",
                                "#f97720", "#b9b9b9", "#c1ad00",
                                "#ad99ff", "#3cd6ff", "#fb85ff",
                                "#f96491"))+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  ggsave('Module_adj_p.pdf',height = 4, width = 3)
  
  ggplot(p_cluster, aes(x = modularity_class, y =BVT_adj_p_pearson, color = modularity_class)) +
    geom_jitter(size = 5)+ theme_bw()+ 
    scale_color_manual(values=c("#2cb384", "#34b5ff", "#38d011",
                                "#f97720", "#b9b9b9", "#c1ad00",
                                "#ad99ff", "#3cd6ff", "#fb85ff",
                                "#f96491"))+
    theme(axis.title = element_text(size = 7), 
          axis.text = element_text(size = 7), 
          legend.text = element_text(size = 7), 
          legend.title = element_text(size = 7))
  ggsave('Module_adj_p_2.pdf',height = 3, width = 3)
  write.csv(p_cluster,'correlated_variables_cluster.csv')
}

# output metadata explanation
{
  setwd('/Users/binzhu/secure/vamp/home/metadata/vahmp/oldDrops')
  data <- read.delim('ClinicalData_Dictionary_full.txt')
  setwd('/Users/binzhu/Desktop/pH/results')
  x = data$DATAFIELD
  for (b in 1: 30) {
    for (a in 1: length(x)) {
      x1= x[a]
      x2 = str_split(x1,'')
      x2 = x2[[1]]
      x2 = x2[length(x2)]
      
      if (!is.na(as.numeric(x2))) {
        x[a] = str_sub(x1, start = 1L, end = -2L)
      }
    }
  }
  data_2 = data
  data_2$DATAFIELD = x
  
  keep = (x %in% metadata_list)
  output = data[keep,]
  write.csv(output,'questions_not_all.csv')
  metadata_list[!metadata_list %in% x]
  
  
  
}

########## application of B-vagitype with DAVT ############
setwd('/Users/binzhu/Desktop/pH/results')
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
  CST_cluster$CST_cluster[CST_cluster$CST_cluster == 9] = 'VI'

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
                               "VI"="#e377c2"),
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

# prepare metadata
{
  
  colnames(metadata_all)[colnames(metadata_all) == 'mom_delivery_method'] = 'mom_delivery_method (C-section yes / vaginal no)'
  colnames(metadata_all)[colnames(metadata_all) == 'Delivery'] = 'Preterm delivery'
  colnames(metadata_all)[colnames(metadata_all) == 'baby_icu'] = 'baby_nicu'
  colnames(metadata_all)[colnames(metadata_all) == 'baby_problems_none'] = 'baby_problems'
  
  metadata_all[metadata_all == 'PTB' | metadata_all == 'Yes_I_plan_to_breast_feed' | 
                 metadata_all == 'Yes_my_baby_was_put_in_the_ICU' |metadata_all == 'Yes_I_have_been_diagnosed_with_hyperemesis' |
                 metadata_all == 'C_Section' |metadata_all == 'Cesarean_section' |  
                 metadata_all == 'Female'| metadata_all == 'Yes_98_100'| metadata_all == 'Yes_101_103'| metadata_all == 'Yes_membranes_ruptured_1_18_hours_before_labor'| metadata_all == 'Yes_membranes_ruptured_before_37_weeks_of_pregnancy'] = 'Yes'
  
  metadata_all[metadata_all == 'TB' | metadata_all == 'No_I_have_not_been_exposed_to_any_chemicals' | metadata_all == 'No_I_have_not_been_to_the_hospital' | metadata_all == 'No_I_have_not_been_diagnosed_with_hyperemesis_and_I_do_not_have_symptoms' |
                 metadata_all == 'No_my_baby_was_not_put_in_the_ICU' | metadata_all == 'No_I_do_not_plan_to_breast_feed' | metadata_all == 'No_I_am_not_taking_any_new_medications' | metadata_all == 'No_I_have_not_been_diagnosed_with_hyperemesis_but_I_have_symptoms' |
                 metadata_all == 'Vaginal' | metadata_all == 'Vaginal_delivery' |  metadata_all == 'No_I_have_not_had_a_fever_since_my_last_visit' |
                 metadata_all == 'Male' | metadata_all == 'none'| metadata_all == 'No_there_has_been_no_premature_rupture_of_membranes'] = 'No'
  
  metadata_all[metadata_all == 'Not_Sure' | 
                 metadata_all == 'Not_sure' | 
                 metadata_all == 'Uncertain' |
                 metadata_all == 'Do_not_know'|
                 metadata_all == 'Nelson 6th Floor'|
                 metadata_all == 'Baby_NICU' | metadata_all == 'NaN' | metadata_all == 'Don_t_know'] = NA
  
  # write.csv(metadata_all,'metadata_all_ori.csv') ### run when output contains sample ID or participant ID
  
  # smoker_second_hand
  metadata_all[metadata_all == 'Never'] = 0
  metadata_all[metadata_all == 'Rarely'] = 1
  metadata_all[metadata_all == 'Almost_every_day'] = 2
  metadata_all[metadata_all == 'Every_day'] = 3
  
  # various worries
  metadata_all[metadata_all == 'Much_Less'] = 1
  metadata_all[metadata_all == 'A_Bit_Less'] = 2
  metadata_all[metadata_all == 'About_Average'] = 3
  metadata_all[metadata_all == 'A_Bit_More'] = 4
  metadata_all[metadata_all == 'Much_More'] = 5
  
  metadata_all[metadata_all == 'Never'] = 0
  metadata_all[metadata_all == 'Almost_Never'] = 1
  metadata_all[metadata_all == 'Sometimes'] = 2
  metadata_all[metadata_all == 'Fairly_Often'] = 3
  metadata_all[metadata_all == 'Very_Often'] = 4
  
  
  # yogurt, milk, cheese, ice_cream
  metadata_all[metadata_all == '1_2_servings_per_day'] = 4
  metadata_all[metadata_all == '3_4_servings_per_day'] = 5
  metadata_all[metadata_all == '5_servings_per_day'] = 6
  metadata_all[metadata_all == 'Less_than_1_serving_per_week'] = 1
  metadata_all[metadata_all == '1_2_servings_per_week'] = 2
  metadata_all[metadata_all == '3_4_servings_per_week'] = 3
  
  # weight_change
  metadata_all[metadata_all == 'My_weight_has_NOT_changed_since_my_last_visit'] = 0
  metadata_all[metadata_all == 'I_have_LOST_weight_since_my_last_visit'] = -1
  metadata_all[metadata_all == 'I_have_GAINED_weight_since_my_last_visit'] = 1
  
  # weight_lost, weight_gained
  metadata_all[metadata_all == '0_5_lbs'] = 1
  metadata_all[metadata_all == '10_15_lbs'] = 3
  metadata_all[metadata_all == '5_10_lbs'] = 2
  metadata_all[metadata_all == 'More_than_15_lbs'] = 4
  
  # physical_activity_vigorous, physical_activity_moderate, physical_activity_light
  metadata_all[metadata_all == '0_times'] = 0
  metadata_all[metadata_all == '1_2_times'] = 1
  metadata_all[metadata_all == '3_4_times'] = 2
  metadata_all[metadata_all == '5_6_times'] = 3
  metadata_all[metadata_all == '7_times'] = 4
  
  # uti_lifetime_number, 
  metadata_all[metadata_all == '2_to_4'] = 3
  metadata_all[metadata_all == '2_to_4'] = 3
  
  # douche_frequency, vaginal_sex_frequency, vaginal_penetration_frequency, received_oral_frequency, gave_oral_frequency, anal_sex_frequency, vaginal_lubrication_frequency
  metadata_all[metadata_all == '1_3_times_per_month'] = 2
  metadata_all[metadata_all == '2_6_times_per_week'] = 4
  metadata_all[metadata_all == 'Less_than_once_a_month'] = 1
  metadata_all[metadata_all == 'Never'] = 0
  metadata_all[metadata_all == 'Once_a_day'] = 5
  metadata_all[metadata_all == 'Once_a_week'] = 3
  
  # income
  metadata_all[metadata_all == '15_000_19_999'] = 2
  metadata_all[metadata_all == '20_000_39_999'] = 3
  metadata_all[metadata_all == '40_000_59_999'] = 4
  metadata_all[metadata_all == '60_000_79_999'] = 5
  metadata_all[metadata_all == '80_000_or_more'] = 6
  metadata_all[metadata_all == 'Under_15_000'] = 1
  
  # education
  metadata_all[metadata_all == '2_year_College_Degree'] = 3
  metadata_all[metadata_all == '4_year_College_Degree'] = 4
  metadata_all[metadata_all == 'Doctoral_or_Professional_Degree'] = 6
  metadata_all[metadata_all == 'High_School_GED'] = 2
  metadata_all[metadata_all == 'Less_than_High_School'] = 1
  metadata_all[metadata_all == 'Masters_Degree'] = 5
  metadata_all[metadata_all == 'Some_College'] = NA
  
  # prenatal_care_start
  metadata_all[metadata_all == '0_4_weeks_after_conception'] = 3
  metadata_all[metadata_all == '4_weeks_to_the_end_of_the_1st_trimester'] = 2
  metadata_all[metadata_all == '4_weeks_to_the_end_of_the_first_trimester'] = 2
  metadata_all[metadata_all == 'Before_conception'] = 4
  metadata_all[metadata_all == 'During_the_2nd_trimester'] = 1
  metadata_all[metadata_all == 'During_the_second_trimester'] = 1
  
  metadata_all[metadata_all == 'Yes'] =1
  metadata_all[metadata_all == 'No'] =0
  
  metadata = metadata_all
  metadata = apply(metadata, 2, as.numeric)
  keep = sapply(1:ncol(metadata), function(j) (length(unique(metadata[,j])) >= 2 ))
  metadata = metadata[,keep]
  
  metadata = as.data.frame(metadata)
  
  for (a in 1: ncol(metadata)) {
    data = metadata[,a]
    keep = unique(data); keep = keep[!is.na(keep)]
    if (length(keep) >2 & length(keep) <=6) {
      keep = data > mean(data[!is.na(data)]) & !is.na(data); metadata[keep,a] = 1
      keep = data <= mean(data[!is.na(data)]) & !is.na(data); metadata[keep,a] = 0
    } else {next}
    
  }
  
  
  metadata_all[metadata_all == 'Yes'] =1
  metadata_all[metadata_all == 'No'] =0
  
  metadata = metadata_all
  metadata = apply(metadata, 2, as.numeric)
  keep = sapply(1:ncol(metadata), function(j) (length(unique(metadata[,j])) >= 2 ))
  metadata = metadata[,keep]
  
  metadata = as.data.frame(metadata)
  
  for (a in 1: ncol(metadata)) {
    data = metadata[,a]
    keep = unique(data); keep = keep[!is.na(keep)]
    if (length(keep) ==2) {
      data[data == 0] = 'No'; data[data == '1'] = 'Yes'; metadata[,a] = data

    } else {next}
    
  }
  
}

# test classification
metadata_test = data.frame(CST = CST, Vagitype= mytypes, BVT = BVT) 
for (n in 1:3) {
  output_table = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = length(unique(metadata_test[,n]))))
  
  colnames(output_table) = colnames(metadata)
  row.names(output_table) = sort(unique(metadata_test[,n]))
  
  output_table_2 = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = 2))
  colnames(output_table_2) = colnames(metadata)
  row.names(output_table_2) = c('Kruskal–Wallis test',"Pearson's Chi-squared test")
  
  output_table_3 = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = 1))
  colnames(output_table_3) = colnames(metadata)
  
  for (y in 1: ncol(output_table)) {
    data = metadata[,y]
    data = data[!is.na(data)]
    if (length(unique(data)) > 6) {
      data_2 = median(data)
      data_2 = round(data_2, 1)
      data_3 = quantile(data)
      data_3 = round(data_3, 1)
      output_table_3[1,y] = paste0(data_2, " (",data_3[2], ", ",data_3[4], "); ",length(data))
    } else if (length(unique(data)) == 2) {
      data_2 = sum(data == 'Yes') / length(data)
      data_2 = round(data_2, 3)
      output_table_3[1,y] = paste0(data_2*100,'%; ',length(data))
    } else {next}
    
    for (x in 1:nrow(output_table)) {
      data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
      data = data[!is.na(data)]
      
      if (length(unique(data)) > 6) {
        data_2 = median(data)
        data_2 = round(data_2, 1)
        data_3 = quantile(data)
        data_3 = round(data_3, 1)
        output_table[x,y] = paste0(data_2, " (",data_3[2], ", ",data_3[4], "); ",length(data))
      } else if (length(unique(data)) == 2) {
        data_2 = sum(data == 'Yes') / length(data)
        data_2 = round(data_2, 3)
        output_table[x,y] = paste0(data_2*100,'%; ',length(data))
      } else {next}
    }
    
    data_4 = as.data.frame(cbind(metadata_test[,n], metadata[,y]))
    keep = !is.na(data_4[,2])
    data_4 = data_4[keep,]
    
    if (length(unique(data_4$V2)) == 2 | length(unique(data_4$V2)) >6) {
      output_table_2[1,y] = kruskal.test(V2~V1, data_4)$p.value
      output_table_2[2,y] = chisq.test(data_4$V2,  data_4$V1, correct=FALSE)$p.value
    }
    
  }
  write.csv(output_table_3,'Sample_size.csv',row.names = F)
  
  output_table = rbind(output_table,output_table_2)

  if (n ==1) {output_table_CST = output_table} else if (n ==2) {output_table_vagitype = output_table}  else {output_table_BVT = output_table} 
}
{
  keep =  colSums(is.na(output_table_CST)) < nrow(output_table_CST)/3 & 
    colSums(is.na(output_table_vagitype)) < nrow(output_table_vagitype)/3 & 
    colSums(is.na(output_table_BVT)) < nrow(output_table_BVT)/3
  
  sum(keep)

  output_table_CST = output_table_CST[ , keep]
  output_table_vagitype = output_table_vagitype[ ,keep]
  output_table_BVT = output_table_BVT[ ,keep]
  
  compare_list = data.frame(Metadata = colnames(output_table_CST), CST_pvalue = as.numeric(as.character(output_table_CST[11,])),
                            Vagitype_pvalue = as.numeric(as.character(output_table_vagitype[15,])), 
                            BVT_pvalue = as.numeric(as.character(output_table_BVT[10,])), Best_pvalue = NA)
  
  for (a in 2:4) {
    data <- adjust.p(compare_list[,a], pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
    data <- data$adjp
    data <- data$adjusted.p
    compare_list[,a] = data
  }
  output_table_CST[11,] = compare_list[,2];output_table_vagitype[15,] = compare_list[,3];output_table_BVT[10,] = compare_list[,4]
  
  compare_list_3 = data.frame(Metadata = colnames(output_table_CST), CST_pvalue = as.numeric(as.character(output_table_CST[10,])),
                            Vagitype_pvalue = as.numeric(as.character(output_table_vagitype[14,])), 
                            BVT_pvalue = as.numeric(as.character(output_table_BVT[9,])), Best_pvalue = NA)
  
  for (a in 2:4) {
    data <- adjust.p(compare_list_3[,a], pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
    data <- data$adjp
    data <- data$adjusted.p
    compare_list_3[,a] = data
  }
  output_table_CST[10,] = compare_list_3[,2];output_table_vagitype[14,] = compare_list_3[,3];output_table_BVT[9,] = compare_list_3[,4]
  

  compare_list_2 = sapply(1:nrow(compare_list), function(j) (compare_list[j,5] = which.min(compare_list[j,c(2:4)])))
  compare_list_2[compare_list_2 == 1] = 'CST';compare_list_2[compare_list_2 == 2] = 'Vagitype';compare_list_2[compare_list_2 == 3] = 'BVT';
  compare_list[,5] = compare_list_2
  
  keep = compare_list[,2] <= 0.05 | compare_list[,3] <= 0.05 | compare_list[,4] <= 0.05; sum(keep)
  compare_list[!keep,5] = ''
  
  metadata_list = compare_list$Metadata[!is.na(compare_list$Best_pvalue)]
#  metadata_list = read.csv('metadata_list.csv', header = F)
#  metadata_list = metadata_list$V1
  output_table_CST = output_table_CST[,colnames(output_table_CST) %in% metadata_list]
  output_table_vagitype = output_table_vagitype[,colnames(output_table_vagitype) %in% metadata_list]
  output_table_BVT = output_table_BVT[,colnames(output_table_BVT) %in% metadata_list]
#  compare_list = compare_list[compare_list$Metadata %in% metadata_list,]
  
  write.csv(output_table_CST, paste0('Distribution_CST.csv'))
  write.csv(output_table_vagitype, paste0('Distribution_Vagitype.csv'))
  write.csv(output_table_BVT, paste0('Distribution_BVT.csv'))
  output_table = rbind(output_table_CST,output_table_vagitype,output_table_BVT)
  write.csv(output_table, paste0('Distribution_all.csv'))
  write.csv(compare_list, paste0('compare_list.csv'))
}

# differential abundance
{
  # select variables
  {
#    da_list = compare_list$Metadata[compare_list$CST_pvalue <= 0.001 & 
#                                      compare_list$Vagitype_pvalue <= 0.001 &
#                                      compare_list$BVT_pvalue <= 0.001]
#    da_list; length(da_list)
#    write.csv(da_list,'da_list.csv',row.names = F)
    
    # find correlated variables
    metadata_list = read.csv('metadata_list.csv',header = F)
    metadata_list = metadata_list$V1
    metadata = metadata[,colnames(metadata) %in% metadata_list]
    data =  t(metadata)
    output = newwork_rcorr(data, normalization_method = NA, type = 'spearman', pvalue = 0.05, 
                           cor_parameter= 0, style = 1, bar_max = 1, bar_min = -1, 
                           pheatmap_fontsize = 5, treeheight = 20, alpha = 0.05)
    p.yes.rr = abs(output$cor_matrix)
    
    bar_max = max(p.yes.rr,na.rm = T)
    bar_min = min(p.yes.rr,na.rm = T)
    paletteLength <- 50
    myColor <- colorRampPalette(c("white", "red"))(paletteLength)
    p = pheatmap(p.yes.rr, color=myColor, fontsize = 5, 
                 treeheight_row = 50, treeheight_col = 50)

    pdf("correlated_variables.pdf", width=12, height=12)
    print(p)
    dev.off()
    
    correlated_variables_network = output$gephi_input
    correlated_variables_network$Weight = abs(correlated_variables_network$Weight)
    correlated_variables_network$Weight = as.numeric(as.character(correlated_variables_network$Weight))
    write.csv(correlated_variables_network,'correlated_variables_network.csv', row.names = F)
    
    data_3 = as.data.frame(output$cor_matrix)
    data_3 = gather(data_3)
    data_4 = as.data.frame(output$adj_p)
    data_4 = gather(data_4)
    data_3$adj_p = data_4$value
    data_3$variable_2 = rep(colnames())
    
    write.csv(output$cor_matrix,'correlated_variables_network_2.csv')
    
    p_cluster = read.csv('gephi_metadata_network.csv')
#    p_cluster = sort(cutree(p$tree_col, k=30))
#    p_cluster = as.data.frame(p_cluster)
    
    p_cluster = p_cluster[order(p_cluster$betweenesscentrality,decreasing = T),]
    table(p_cluster$modularity_class)
    modul_list = unique(p_cluster$modularity_class)
    
    for (a in 1: length(modul_list)) {
      n = which(p_cluster$modularity_class == modul_list[a])
      p_cluster$modularity_class[n[4:length(n)]] = NA
    }
    p_cluster = p_cluster[!is.na(p_cluster$modularity_class),]
    
    write.csv(p_cluster,'correlated_variables_cluster.csv')
    
    da_list = p_cluster$Label
  }

  # differencial abundance analysis
  {
    for (a in 1:length(da_list)) {
      if (da_list[a] %in% data_all$Metadata) {next}
      
      metadata_2 = metadata[,which(colnames(metadata)==da_list[a])]
      keep = !is.na(metadata_2)
      metadata_2 = metadata_2[keep]
      reads_table_2 = reads_table_all[,keep]
      
      if (length(unique(metadata_2)) >2){
        keep = metadata_2 > median(metadata_2)
        metadata_2[keep] = 1
        metadata_2[!keep] = 0
      }
      
      data = dif_abundance2(reads_table_2,metadata_2)
      data = data$data
      
      data$Metadata = da_list[a]
      data$Taxa = row.names(data)
      
      if (a ==1) {
        data_all = data 
      } else {data_all = rbind(data_all,data)}
    }
    
    write.csv(data_all,'differential_abundance.csv')
    
    data = data_all[data_all$we.eBH <= 0.05,]
    data$log_adj_pvalue = -log10(data$we.eBH)
    ggplot(data, aes(Taxa, Metadata)) + 
      geom_point(aes(col= diff.btw, size=log_adj_pvalue)) + 
      scale_color_gradient2(midpoint=0, low="blue", mid="white",
                            high="red", space ="Lab" ) +
      coord_flip() +theme_bw()+          # convert x y axis
      labs(x = 'Taxa')+ 
      theme(axis.title = element_text(size = 12), 
            axis.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 12))
    ggsave('differential_abundance.pdf',width=8, height=8)
  }
  
  # New classification DAVT
  {
    data = data[,c(8,12,13)]
    
    metadata_list_2 = unique(data$Metadata); taxa_list_2 = unique(data$Taxa)
    temp = as.data.frame(matrix(data = 0, nrow = length(taxa_list_2), ncol = length(metadata_list_2))) 
    colnames(temp) = metadata_list_2; row.names(temp) = taxa_list_2
    
    for (a in 1: nrow(data)) {
      x = which(data$Taxa[a] == taxa_list_2); y = which(data$Metadata[a] == metadata_list_2)
      temp[x,y]= data$diff.btw[a]
    }
    
    taxa_list_2 <- as.data.frame(taxa_list_2)
    row.names(taxa_list_2) = row.names(temp)
    taxa_cluster = read.csv('/Users/binzhu/Desktop/pH/Cluster.csv')
    taxa_list_2$BVT = NA
    
    for (a in 1:nrow(taxa_list_2)) {
      n= which(taxa_cluster$Taxa ==taxa_list_2$taxa_list_2[a])
      if (length(n) ==1) {
        taxa_list_2$BVT[a] = taxa_cluster$Cluster[n]
      } else {taxa_list_2$BVT[a] = 'VIII'}
      
    }
    taxa_list_2$taxa_list_2 = NULL

    bar_max = max(temp,na.rm = T)
    bar_min = min(temp,na.rm = T)
    paletteLength <- 50
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
    myBreaks <- c(seq(bar_min, 0, length.out=ceiling(paletteLength/2) +1), 
                  seq(max(temp,na.rm = T)/paletteLength, bar_max, length.out=floor(paletteLength/2)))

    p = pheatmap(temp,cluster_rows = T, color=myColor, cluster_cols =T, clustering_distance_cols = "minkowski", clustering_method = "ward.D",
                 treeheight_col = 50)

    p_cluster_2 = sort(cutree(p$tree_row, k=6))
    p_cluster_2 = as.data.frame(p_cluster_2)
    p_cluster_2$Taxa = row.names(p_cluster_2)
    p_cluster_2 = p_cluster_2[match(row.names(taxa_list_2), row.names(p_cluster_2)),]
    taxa_list_2$DAVT = p_cluster_2$p_cluster_2
    taxa_list_2$DAVT[taxa_list_2$DAVT == 1] = 'I';taxa_list_2$DAVT[taxa_list_2$DAVT ==2] = 'II'
    taxa_list_2$DAVT[taxa_list_2$DAVT == 3] = 'VI';taxa_list_2$DAVT[taxa_list_2$DAVT ==4] = 'V'
    taxa_list_2$DAVT[taxa_list_2$DAVT == 5] = 'IV';taxa_list_2$DAVT[taxa_list_2$DAVT ==6] = 'III'
    #taxa_list_2$DAVT[taxa_list_2$DAVT == 7] = 'V';taxa_list_2$DAVT[taxa_list_2$DAVT ==8] = 'IV'
    
    pdf("differential_abundance_heatmap.pdf", width=7, height=8)
    print(pheatmap(temp,cluster_rows = T, color=myColor, breaks=myBreaks, cluster_cols =T, 
                   clustering_distance_cols = "minkowski", clustering_method = "ward.D",
                   annotation_row = taxa_list_2, treeheight_col = 50))
    dev.off()
    
    mycolors = list(BVT = c("I"="#FFDC00","II"="#49C800","III"="#aec7e8","IV"="#17becf","V"="#7f7f7f","VI"="#ff7f0e","VII"="#d62728","VIII"="#9C9146"),
                    DAVT = c("I"="#FFDC00","II"="#c5b0d5","III"="#aec7e8","IV"="#17becf","V"="#d62728","VI"="#ff7f0e"))
    
    taxa_list_2$BVT = factor(taxa_list_2$BVT, levels = c('I','II',"III",'IV','V','VI','VII','VIII'))
    taxa_list_2$DAVT = factor(taxa_list_2$DAVT, levels = c('I','II',"III",'IV','V','VI'))

    library(grid)
    draw_colnames_45 <- function (coln, gaps, ...) {
      coord = pheatmap:::find_coordinates(length(coln), gaps)
      x = coord$coord - 0.5 * coord$size
      res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
      return(res)}
    
    assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                      ns=asNamespace("pheatmap"))
    
    pdf("differential_abundance_heatmap.pdf", width=7, height=8)
    print(pheatmap(temp,cluster_rows = T, color=myColor, breaks=myBreaks, cluster_cols =T, 
                   clustering_distance_cols = "minkowski", clustering_method = "ward.D",
                   annotation_row = taxa_list_2, annotation_colors =mycolors,treeheight_col = 50))
    dev.off()
    
    library(grid)
    draw_colnames_45 <- function (coln, gaps, ...) {
      coord = pheatmap:::find_coordinates(length(coln), gaps)
      x = coord$coord - 0.5 * coord$size
      res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
      return(res)}
    
    assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                      ns=asNamespace("pheatmap"))
    
    pdf("differential_abundance_heatmap_2.pdf", width=7, height=8)
    print(pheatmap(temp,cluster_rows = T, color=myColor, breaks=myBreaks, cluster_cols =T, 
                   clustering_distance_cols = "minkowski", clustering_method = "ward.D",treeheight_col = 50))
    dev.off()
    
    write.csv(taxa_list_2,'DAVT_cluster.csv')
  }

  # get DAVT
  # get biological vagitype (B-vagitype)
  {
    # get DAVT
    {
      type_th = 0.3
      
      reads_table = as.data.frame(t(reads_table_all_abundance))
      reads_table <- reads_table[,apply(reads_table,2,sum) > 0]
      mytypes <- apply(reads_table,1,which.max)   # find the most abundant taxa
      maxprop <- reads_table[matrix(c(1:nrow(reads_table),mytypes), ncol=2)]  # the abundance of the most abundant taxa
      mytypes <- colnames(reads_table)[mytypes]   # find the name of the most abundant taxa
      mytypes[maxprop < type_th] <- "No_Type"
      
      DAVT = as.data.frame(mytypes)
      taxa_order = taxa_list_2
      taxa_order$BVT = NULL
      taxa_order$Taxa = row.names(taxa_order)
      taxa_order$DAVT = as.character(taxa_order$DAVT)
      taxa_order = rbind(taxa_order, c('Others','No_Type'))
      colnames(taxa_order)[1] = 'Cluster'
      
      n = match(DAVT$mytypes, taxa_order$Taxa)
      DAVT$DAVT = taxa_order$Cluster[n]
      DAVT = DAVT$DAVT
    }
    
    # get vagitype
    {
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
      
      keep = mytypes %in% top_vagitype$mytypes[1:12]
      mytypes[!keep] = 'Others'
    }
  }
  
  # test classification
  {
    metadata_test = data.frame(CST = CST, Vagitype= mytypes, BVT = BVT, DAVT = DAVT) 
    write.csv(metadata_test,'VMB_classification.csv')
    for (n in 1:4) {
      output_table = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = length(unique(metadata_test[,n]))))
      
      colnames(output_table) = colnames(metadata)
      row.names(output_table) = sort(unique(metadata_test[,n]))
      
      output_table_2 = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = 2))
      colnames(output_table_2) = colnames(metadata)
      row.names(output_table_2) = c('Kruskal–Wallis test',"Pearson's Chi-squared test")
      
      
      for (y in 1: ncol(output_table)) {
        for (x in 1:nrow(output_table)) {
          data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
          data = data[!is.na(data)]
          
          if (length(unique(data)) > 6) {
            data_2 = median(data)
            data_2 = round(data_2, 1)
            data_3 = quantile(data)
            data_3 = round(data_3, 1)
            output_table[x,y] = paste0(data_2, " (",data_3[2], ", ",data_3[4], "); ",length(data))
          } else if (length(unique(data)) == 2) {
            data_2 = sum(data == 1) / length(data)
            data_2 = round(data_2, 3)
            output_table[x,y] = paste0(data_2*100,'%; ',length(data))
          } else {next}
        }
        
        data_4 = as.data.frame(cbind(metadata_test[,n], metadata[,y]))
        keep = !is.na(data_4[,2])
        data_4 = data_4[keep,]
        
        if (length(unique(data_4$V2)) == 2 | length(unique(data_4$V2)) >6) {
          output_table_2[1,y] = kruskal.test(V2~V1, data_4)$p.value
          output_table_2[2,y] = chisq.test(data_4$V2,  data_4$V1, correct=FALSE)$p.value
        }
        
      }
      output_table = rbind(output_table,output_table_2)
      
      if (n ==1) {output_table_CST = output_table} else if (n ==2) 
      {output_table_vagitype = output_table}  else if (n == 3) 
          {output_table_BVT = output_table} else (output_table_DAVT = output_table)
    }
    {
      keep =  colSums(is.na(output_table_CST)) < nrow(output_table_CST)/2 & 
        colSums(is.na(output_table_vagitype)) < nrow(output_table_vagitype)/2 & 
        colSums(is.na(output_table_BVT)) < nrow(output_table_BVT)/2& 
        colSums(is.na(output_table_DAVT)) < nrow(output_table_DAVT)/2
      
      sum(keep)
      
      output_table_CST = output_table_CST[ , keep]
      output_table_vagitype = output_table_vagitype[ ,keep]
      output_table_BVT = output_table_BVT[ ,keep]
      output_table_DAVT = output_table_DAVT[ ,keep]
      
      compare_list = data.frame(Metadata = colnames(output_table_CST), CST_pvalue = as.numeric(as.character(output_table_CST[11,])),
                                Vagitype_pvalue = as.numeric(as.character(output_table_vagitype[15,])), 
                                BVT_pvalue = as.numeric(as.character(output_table_BVT[10,])), 
                                DAVT_pvalue = as.numeric(as.character(output_table_DAVT[9,])),Best_pvalue = NA)
      
      for (a in 2:5) {
        data <- adjust.p(compare_list[,a], pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
        data <- data$adjp
        data <- data$adjusted.p
        compare_list[,a] = data
      }
      output_table_CST[11,] = compare_list[,2];output_table_vagitype[15,] = compare_list[,3];
      output_table_BVT[10,] = compare_list[,4];output_table_DAVT[9,] = compare_list[,5];
      
      compare_list_3 = data.frame(Metadata = colnames(output_table_CST), CST_pvalue = as.numeric(as.character(output_table_CST[10,])),
                                  Vagitype_pvalue = as.numeric(as.character(output_table_vagitype[14,])), 
                                  BVT_pvalue = as.numeric(as.character(output_table_BVT[9,])), 
                                  DAVT_pvalue = as.numeric(as.character(output_table_DAVT[8,])), Best_pvalue = NA)
      
      for (a in 2:5) {
        data <- adjust.p(compare_list_3[,a], pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
        data <- data$adjp
        data <- data$adjusted.p
        compare_list_3[,a] = data
      }
      output_table_CST[10,] = compare_list_3[,2];output_table_vagitype[14,] = compare_list_3[,3];
      output_table_BVT[9,] = compare_list_3[,4];output_table_DAVT[8,] = compare_list_3[,5];
      
      
      compare_list_2 = sapply(1:nrow(compare_list), function(j) (compare_list[j,6] = which.min(compare_list[j,c(2:5)])))
      compare_list_2[compare_list_2 == 1] = 'CST';compare_list_2[compare_list_2 == 2] = 'Vagitype';
      compare_list_2[compare_list_2 == 3] = 'BVT';compare_list_2[compare_list_2 == 4] = 'DAVT';
      compare_list[,6] = compare_list_2
      
      keep = compare_list[,2] <= 0.05 | compare_list[,3] <= 0.05 | compare_list[,4] <= 0.05 | compare_list[,5] <= 0.05; sum(keep)
      compare_list[!keep,6] = ''
      
      metadata_list = compare_list$Metadata[!is.na(compare_list$Best_pvalue)]
      #  metadata_list = read.csv('metadata_list.csv', header = F)
      #  metadata_list = metadata_list$V1
      output_table_CST = output_table_CST[,colnames(output_table_CST) %in% metadata_list]
      output_table_vagitype = output_table_vagitype[,colnames(output_table_vagitype) %in% metadata_list]
      output_table_BVT = output_table_BVT[,colnames(output_table_BVT) %in% metadata_list]
      output_table_DAVT = output_table_DAVT[,colnames(output_table_DAVT) %in% metadata_list]
      #  compare_list = compare_list[compare_list$Metadata %in% metadata_list,]
      
      output_table = rbind(output_table_CST,output_table_vagitype,output_table_BVT,output_table_DAVT)
      write.csv(output_table, paste0('Distribution_all.csv'))
      write.csv(compare_list, paste0('compare_list.csv'))
      
      data = compare_list[,c(2:5)]
      data = gather(data)
      data$key = str_remove_all(data$key,'_pvalue')
      colnames(data) = c('Method','-Log10(Adj-P-value)')
      data$`-Log10(Adj-P-value)` = -log10(data$`-Log10(Adj-P-value)`)

      data$Method = factor(data$Method, levels = c('CST', 'Vagitype', 'BVT', 'DAVT'))
      ggplot(data, aes(x = Method, y = `-Log10(Adj-P-value)`, color = Method)) + 
        scale_y_continuous(breaks = seq(0, 110, 5))+
        geom_boxplot(outlier.shape = NA)+
        geom_jitter(size = 0.1)+ theme_bw()+
        theme(axis.title = element_text(size = 7), 
              axis.text = element_text(size = 7), 
              legend.text = element_text(size = 7), 
              legend.title = element_text(size = 7))
      ggsave('Method.compare.pdf',height = 1, width = 3.5)
      ggsave('Method.compare_2.pdf',height = 8.5, width = 3.5)
      
      # calculate significance
      {
        data$Method = as.character(data$Method)
        factor_levels = unique(data$Method)
        n = length(factor_levels)
        
        Method_sig = as.data.frame(matrix(data = NA, nrow =n, ncol = n))
        colnames(Method_sig) = factor_levels
        row.names(Method_sig) = factor_levels

        for (a in 1:(n-1)) {
          for (b in (a+1) : n) {
            factor_level1 <- data$`-Log10(Adj-P-value)`[data$Method == factor_levels[a]]
            factor_level2 <- data$`-Log10(Adj-P-value)`[data$Method == factor_levels[b]]
            
            Method_sig[a,b] <- wilcox.test(factor_level1, factor_level2, paired = T)$p.value
            
          }
        }
        write.csv(Method_sig,'Method.compare.csv')
      }

    }
  }
  
  # heatmap
  {
    Dominant_species_2$DAVT = DAVT
    Dominant_species_2$Vagitype = str_replace_all(Dominant_species_2$Vagitype, '_',' ')
    mycolors = list(Vagitype = c("Lactobacillus crispatus"="#FFDC00","Lactobacillus gasseri"="#c5b0d5","Lactobacillus iners"="#aec7e8",
                                 "Lactobacillus jensenii"="#49C800","Ca. L. vaginae"="#ff7f0e","Gardnerella vaginalis"="#d62728","No Type"="#4d4d4d",
                                 "Atopobium vaginae"="#c49c94","Prevotella bivia"="#17becf","Mycoplasma"="#e377c2",
                                 "Sneathia amnii"="#FF00D1", "Streptococcus"="#7f7f7f","Others"="#9C9146"),
                    CST = c("I"="#FFDC00","II"="#c5b0d5","III"="#aec7e8","V"="#49C800","IVA"="#ff7f0e","IVB"="#d62728","IVC"="#4d4d4d","IVD"="#17becf",
                            "VI"="#e377c2"),
                    BVT = c("I"="#FFDC00","II"="#49C800","III"="#aec7e8","VI"="#ff7f0e","VII"="#d62728","IV"="#17becf","V"="#7f7f7f",
                            "Others"="#9C9146"),
                    DAVT = c("I"="#FFDC00","II"="#c5b0d5","III"="#aec7e8","IV"="#17becf","V"="#d62728","VI"="#ff7f0e","Others"="#9C9146"))
    
    Dominant_species_2 = Dominant_species_2[,c(4,1,2,3)]
    pdf("heatmap_clustering.pdf", width=18, height=12)
    print(pheatmap(reads_table_cluster,cluster_rows = T, cluster_cols =T, clustering_distance_cols = "euclidean", clustering_method = "complete",
                   annotation_col = Dominant_species_2, annotation_colors =mycolors,show_colnames = F,treeheight_col = 50))
    dev.off()
    
    pdf("heatmap_clustering_legend.pdf", width=18, height=12)
    print(pheatmap(reads_table_cluster,cluster_rows = T, cluster_cols =T, clustering_distance_cols = "euclidean", clustering_method = "complete",
                   annotation_col = Dominant_species_2, annotation_colors =mycolors,show_colnames = F,treeheight_col = 50))
    dev.off()
  }
}

########## application of B-vagitype 5 types ############
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
  
  data = Dominant_species_2[,-3]
  pdf("heatmap_clustering_5_CST_5_vagitype_2.pdf", width=18, height=8)
  print(pheatmap(reads_table_cluster,cluster_rows = T, cluster_cols =T, clustering_distance_cols = "euclidean", clustering_method = "complete",
                 annotation_col = data,show_colnames = F,treeheight_col = 50, annotation_colors =mycolors))
  dev.off()
}

# compare second version no CST_2
metadata_test = Dominant_species_2[,c(4,2,1)]; colnames(metadata_test)[1] = 'CST'
for (n in 1:3) {
  output_table = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = length(unique(metadata_test[,n]))))
  
  colnames(output_table) = colnames(metadata)
  row.names(output_table) = sort(unique(metadata_test[,n]))
  
  output_table_2 = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = 2))
  colnames(output_table_2) = colnames(metadata)
  row.names(output_table_2) = c('Kruskal–Wallis test',"Pearson's Chi-squared test")
  
  
  for (y in 1: ncol(output_table)) {
    for (x in 1:nrow(output_table)) {
      data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
      data = data[!is.na(data)]
      
      if (length(unique(data)) > 6) {
        data_2 = median(data)
        data_2 = round(data_2, 1)
        data_3 = quantile(data)
        data_3 = round(data_3, 1)
        output_table[x,y] = paste0(data_2, " (",data_3[2], ", ",data_3[4], "); ",length(data))
      } else if (length(unique(data)) == 2) {
        data_2 = sum(data == 1) / length(data)
        data_2 = round(data_2, 3)
        output_table[x,y] = paste0(data_2*100,'%; ',length(data))
      } else {next}
    }
    
    data_4 = as.data.frame(cbind(metadata_test[,n], metadata[,y]))
    keep = !is.na(data_4[,2])
    data_4 = data_4[keep,]
    
    if (length(unique(data_4$V2)) == 2 | length(unique(data_4$V2)) >6) {
      output_table_2[1,y] = kruskal.test(V2~V1, data_4)$p.value
      output_table_2[2,y] = chisq.test(data_4$V2,  data_4$V1, correct=FALSE)$p.value
    }
    
  }
  output_table = rbind(output_table,output_table_2)
  
  if (n ==1) {output_table_CST = output_table} else if (n ==2) {output_table_vagitype = output_table}  else {output_table_BVT = output_table} 
}
{
  keep =  colSums(is.na(output_table_CST)) < nrow(output_table_CST)/3 & 
    colSums(is.na(output_table_vagitype)) < nrow(output_table_vagitype)/3 & 
    colSums(is.na(output_table_BVT)) < nrow(output_table_BVT)/3
  
  sum(keep)
  
  output_table_CST = output_table_CST[ , keep]
  output_table_vagitype = output_table_vagitype[ ,keep]
  output_table_BVT = output_table_BVT[ ,keep]
  
  compare_list = data.frame(Metadata = colnames(output_table_CST), CST_pvalue = as.numeric(as.character(output_table_CST[7,])),
                            Vagitype_pvalue = as.numeric(as.character(output_table_vagitype[7,])), 
                            BVT_pvalue = as.numeric(as.character(output_table_BVT[7,])), Best_pvalue = NA)
  
  for (a in 2:4) {
    data <- adjust.p(compare_list[,a], pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
    data <- data$adjp
    data <- data$adjusted.p
    compare_list[,a] = data
  }
  
  output_table_CST[7,] = compare_list[,2];output_table_vagitype[7,] = compare_list[,3];output_table_BVT[7,] = compare_list[,4]
  
  compare_list_3 = data.frame(Metadata = colnames(output_table_CST), CST_pvalue = as.numeric(as.character(output_table_CST[6,])),
                              Vagitype_pvalue = as.numeric(as.character(output_table_vagitype[6,])), 
                              BVT_pvalue = as.numeric(as.character(output_table_BVT[6,])), Best_pvalue = NA)
  
  for (a in 2:4) {
    data <- adjust.p(compare_list_3[,a], pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
    data <- data$adjp
    data <- data$adjusted.p
    compare_list_3[,a] = data
  }
  output_table_CST[6,] = compare_list_3[,2];output_table_vagitype[6,] = compare_list_3[,3];output_table_BVT[6,] = compare_list_3[,4]
  
  
  compare_list_2 = sapply(1:nrow(compare_list), function(j) (compare_list[j,5] = which.min(compare_list[j,c(2:4)])))
  compare_list_2[compare_list_2 == 1] = 'CST';compare_list_2[compare_list_2 == 2] = 'Vagitype';compare_list_2[compare_list_2 == 3] = 'BVT';
  compare_list[,5] = compare_list_2
  
  keep = compare_list[,2] <= 0.05 | compare_list[,3] <= 0.05 | compare_list[,4] <= 0.05; sum(keep)
  compare_list[!keep,5] = NA
  
  metadata_list = compare_list$Metadata[!is.na(compare_list$Best_pvalue)]
  #  metadata_list = read.csv('metadata_list.csv', header = F)
  #  metadata_list = metadata_list$V1
  output_table_CST = output_table_CST[,colnames(output_table_CST) %in% metadata_list]
  output_table_vagitype = output_table_vagitype[,colnames(output_table_vagitype) %in% metadata_list]
  output_table_BVT = output_table_BVT[,colnames(output_table_BVT) %in% metadata_list]
  #  compare_list = compare_list[compare_list$Metadata %in% metadata_list,]
  
  write.csv(output_table_CST, paste0('Distribution_CST_5_types.csv'))
  write.csv(output_table_vagitype, paste0('Distribution_Vagitype_5_types.csv'))
  write.csv(output_table_BVT, paste0('Distribution_BVT_5_types.csv'))
  output_table = rbind(output_table_CST,output_table_vagitype,output_table_BVT)
  write.csv(output_table, paste0('Distribution_all_5_types.csv'))
  write.csv(compare_list, paste0('compare_list_5_types.csv'))
}


# compare first version with CST_2 not used
{
  metadata_test = Dominant_species_2
  for (n in 1:4) {
    output_table = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = 5))
    colnames(output_table) = colnames(metadata)
    row.names(output_table) = sort(unique(metadata_test[,n]))
    
    output_table_2 = as.data.frame(matrix(data = NA, ncol = ncol(metadata), nrow = 2))
    colnames(output_table_2) = colnames(metadata)
    row.names(output_table_2) = c('Kruskal_Wallis test',"Pearson's Chi-squared test")
    
    for (y in 1: ncol(output_table)) {
      for (x in 1:nrow(output_table)) {
        data = metadata[metadata_test[,n] == row.names(output_table)[x],y]
        data = data[!is.na(data)]
        if (length(unique(data)) > 2) {
          data_2 = median(data)
          data_2 = round(data_2, 1)
          data_3 = quantile(data)
          data_3 = round(data_3, 1)
          output_table[x,y] = paste0(data_2, " (",data_3[2], ", ",data_3[4], "); ",length(data))
        } else if (length(unique(data)) == 2) {
          data_2 = sum(data == 1) / length(data)
          data_2 = round(data_2, 3)
          output_table[x,y] = paste0(data_2*100,'%; ',length(data))
        } else {next}
        
      }
      
      data_4 = as.data.frame(cbind(as.character(metadata_test[,n]), metadata[,y]))
      keep = !is.na(data_4[,2])
      data_4 = data_4[keep,]
      
      if (length(unique(data_4$V2)) == 2 | length(unique(data_4$V2)) >2) {
        output_table_2[1,y] = kruskal.test(V2~V1, data_4)$p.value
        output_table_2[2,y] = chisq.test(data_4$V2,  data_4$V1, correct=FALSE)$p.value
      }
      
    }
    output_table = rbind(output_table,output_table_2)
    if (n ==1) {output_table_BVT = output_table} else if (n ==2) {output_table_vagitype = output_table} else if (n ==3) {output_table_CST_2 = output_table} else {output_table_CST_1 = output_table}
  }
  {
    
    keep =  colSums(is.na(output_table_CST_1)) < nrow(output_table_CST_1)/2 & 
      colSums(is.na(output_table_vagitype)) < nrow(output_table_vagitype)/2 & 
      colSums(is.na(output_table_BVT)) < nrow(output_table_BVT)/2 &
      colSums(is.na(output_table_CST_2)) < nrow(output_table_CST_2)/2
    sum(keep)
    
    output_table_CST_1 = output_table_CST_1[ , keep]
    output_table_CST_2 = output_table_CST_2[ , keep]
    output_table_vagitype = output_table_vagitype[ ,keep]
    output_table_BVT = output_table_BVT[ ,keep]
    
    compare_list = data.frame(Metadata = colnames(output_table_CST_1), CST_1_pvalue = as.numeric(as.character(output_table_CST_1[7,])),
                              CST_2_pvalue = as.numeric(as.character(output_table_CST_2[7,])),
                              Vagitype_pvalue = as.numeric(as.character(output_table_vagitype[7,])), 
                              BVT_pvalue = as.numeric(as.character(output_table_BVT[7,])), Best_pvalue = NA)
    
    for (a in 2:5) {
      data <- adjust.p(compare_list[,a], pi0.method="bky", alpha = 0.05,pz = 0.05)  # set alpha for fdr. Another function - p.adjust - use alpha = 0.25
      data <- data$adjp
      data <- data$adjusted.p
      compare_list[,a] = data
    }
    
    
    keep = compare_list[,2] <= 0.05 | compare_list[,3] <= 0.05 | compare_list[,4] <= 0.05 | compare_list[,5] <= 0.05; sum(keep)
    compare_list = compare_list[keep,]; output_table_CST_1 = output_table_CST_1[,keep];output_table_CST_2 = output_table_CST_2[,keep];
    output_table_vagitype = output_table_vagitype[,keep]; output_table_BVT = output_table_BVT[,keep]
    
    
    compare_list_2 = sapply(1:nrow(compare_list), function(j) (compare_list[j,6] = which.min(compare_list[j,c(2:5)])))
    compare_list_2[compare_list_2 == 1] = 'CST_1';compare_list_2[compare_list_2 == 2] = 'CST_2'; compare_list_2[compare_list_2 == 3] = 'Vagitype';compare_list_2[compare_list_2 == 4] = 'BVT';
    compare_list[,6] = compare_list_2
    
    metadata_list = read.csv('metadata_list.csv', header = F)
    metadata_list = metadata_list$V1
    output_table_CST_1 = output_table_CST_1[,colnames(output_table_CST_1) %in% metadata_list]
    output_table_CST_2 = output_table_CST_2[,colnames(output_table_CST_2) %in% metadata_list]
    output_table_vagitype = output_table_vagitype[,colnames(output_table_vagitype) %in% metadata_list]
    output_table_BVT = output_table_BVT[,colnames(output_table_BVT) %in% metadata_list]
    compare_list = compare_list[compare_list$Metadata %in% metadata_list,]
    
    write.csv(output_table_CST_1, paste0('Distribution_CST_1_5type.csv'))
    write.csv(output_table_CST_2, paste0('Distribution_CST_2_5type.csv'))
    write.csv(output_table_vagitype, paste0('Distribution_Vagitype_5type.csv'))
    write.csv(output_table_BVT, paste0('Distribution_BVT_5type.csv'))
    write.csv(compare_list, paste0('compare_list_5type.csv'))
  }
}

