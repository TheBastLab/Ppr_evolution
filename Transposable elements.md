# FasTE approach  
###for help see https://github.com/ellenbell/FasTE  

## Step1 creating TE library  
create a species specific TE library with EDTA 

## Step2 refine classifications with DeepTE  
### prerequisites  
make conda env 'DeepTE' with python version 3.7.1, installed:  
biopython 1.78  
numpy 1.16.0  
tensorflow 1.14.0  
keras 2.2.4  
Download the latest DeepTE scripts.  
### Use DeepTE to classify transposons in the libraries into families
```sh
mkdir working_dir  
mkdir output_dir  
```
```sh
python /NVME/Software/DeepTE/DeepTE.py -d working_dir -o output_dir -i /path/to/TElib.fa -m_dir /path/to/Metazoans_model/ -sp M -fam All
```
  -d Working directory to store intermediate files of each step.  
  -o Output directory to store the output files.  
  -i Input sequences that are unknown TE or DNA sequences.  
  -m (change to -m_dir) Provide model_dir that could be downloaded from website (optional requirements).  
                        Here Metazoans_model (remember to unpack)  
                        If users set -UNS yes, please provide UNS_model directory that can be downlowed in the above link   
  -sp P or M or F or O. P:Plants, M:Metazoans, F:Fungi, and O: Others  
  -fam Provide TE family name for the input te sequence  
  
## Step3 Clean-up, Header simplification
###taken from ellenbell/FasTE without adaption, worked just fine
```sh
sed -e 's/\(#\).*\(__\)/\1\2/'  [path to DeepTE.fasta] > [path to cleaned up library]
```
  
## Step4 RepeatMasker
in conda env repeatM (RepeatMasker version 4.1.2.p1)
```sh
mkdir repeatmasker_opt_dir
```
essential command
```sh
RepeatMasker /path/to/genome_assembly.fa -lib /path/to/opt_DeepTE.cleaned.fasta -dir repeatmasker_opt_dir
```
###custom options, (including recommendations from FasTE pipeline)  
  -pa(rallel)  
	number of threads being used, give 30 - 40 depending on MOTOKO's workload  

  -a(lignments)  
	Writes alignments in .align output file  
 
  -s  
	Slow search; 0-5% more sensitive, 2-3 times slower than default  

  -nolow  
	Does not mask low_complexity DNA or simple repeats  

  -no_is  
	Skips bacterial insertion element check  


```sh
RepeatMasker -pa 35 -a -s -nolow -no_is -dir repeatmasker_opt_dir -lib /path/to/opt_DeepTE.cleaned.fasta /path/to/Genome_assembly_hap0.fasta
```
  
## Step5 RepeatMasker Clean-up
###taken from ellenbell/FasTE without adaption
```sh
awk '!/\*/' [repeatmasker.out] > [noasterisk_repeatmasker.out]
```
and
```sh
sed 's/-int//' [noasterisk_repeatmasker.out] > [tidy_noasterisk_repeatmasker.out]
```

## Step6 Refinement with RM_Trips
###with R using R version 4.1.2  
###just change the inputs i,j,k,l (it helps to install some dependencies seperately)  
###here is the copied RM_Trips script, better get it from https://github.com/ellenbell/FasTE/blob/main/RM_TRIPS.R  


### RM_TRIPS
#Author = Chris Butler  
#Email = c.butler@uea.ac.uk  
#Date = 22/08/19  


rm(list = ls())  

###set up inputs  
i <- '/home/jens/Viki/stuff_for_RM_Trips' #directory where .out file is located  
j <- "tidy_noasterisk_Hls.hap0.fasta.out" #set name of file  
k <- '/home/jens/Viki/stuff_for_RM_Trips' #directory where the repeatmasker library is found (.lib/fasta file)  
l <- 'opt_DeepTE.cleaned.fasta' #set name of .lib file  

#LOAD PACKAGES - ensure these are installed beforehand  
library(dplyr)  
library(ggplot2)  
library(biomartr)  
library(seqinr)  


#SET UP FUNCTIONS  
read_rm1 <- function(file) {
  rm_file <- readr::read_lines(file = file, skip = 2) ##### changed from skip = 3 as this removed the first entry
  rm_file <- lapply(rm_file, function(x) {
    str.res <- unlist(stringr::str_split(x, "\\s+"))[-1]
    str.res <- str.res[1:14]
    return(str.res)
  })
  rm_file <- tibble::as_tibble(do.call(rbind, rm_file))  
  colnames(rm_file) <- c(
    "sw_score",
    "perc_div",
    "perc_del",
    "perc_insert",
    "qry_id",
    "qry_start",
    "qry_end",
    "qry_left",
    "matching_repeat",
    "repeat_id",
    "matching_class",
    "no_bp_in_complement",
    "in_repeat_start",
    "in_repeat_end"
  )
  qry_end <- qry_start <- NULL  
  
  nrow_before_filtering <- nrow(rm_file)  
  
  suppressWarnings(rm_file <- dplyr::mutate(rm_file,
                                            qry_start = as.integer(qry_start),
                                            qry_end = as.integer(qry_end)))  
  
  
  rm_file <-
    dplyr::filter(
      rm_file,
      !is.na(qry_start),
      !is.na(qry_end)
    )  
  
  rm_file <-
    dplyr::mutate(
      rm_file,
      qry_width = as.integer(qry_end - qry_start + 1L))  
  
  nrow_after_filtering <- nrow(rm_file)  
  
  
  if ((nrow_before_filtering - nrow_after_filtering) > 0)  
    message((nrow_before_filtering - nrow_after_filtering) + 1 , " out of ",nrow_before_filtering," rows ~ ", round(((nrow_before_filtering - nrow_after_filtering) + 1) / nrow_before_filtering, 3) , "% were removed from the imported RepeatMasker file, ",
            "because they contained 'NA' values in either 'qry_start' or 'qry_end'.")  
  
  return(rm_file)
}  



###load input files  
setwd(i)  
RM <- read_rm1(j)  

setwd(k)  
TEsequence <- read.fasta(file = l, seqtype = c("DNA"))  


##########  
#FILTERING STEPS  
##########  

#1) Remove simple repeats  
noisonosim <- RM  
noisonosim <- noisonosim %>% filter(matching_class != "Simple_repeat") #no simple repeats  
noisonosim <- noisonosim %>% filter(matching_class != "Low_complexity") #no low complexity repeats  
noisonosim <- noisonosim %>% filter(matching_class != "Satellite") #no satellites  
noisonosim <- noisonosim %>% filter(matching_class != "rRNA") #no rRNA  
noisonosim <- noisonosim %>% filter(matching_class != "snRNA") #no snRNA  
noisonosim <- noisonosim %>% filter(matching_class != "tRNA") #no tRNA  
noisonosim <- noisonosim %>% filter(matching_class != "ARTEFACT") #no artefacts  





#2) Merge fragments   
#Extract the length of every  reference sequence and store the results in a vector.  
seq_length <- rep(NA,length(TEsequence))  
for (f in 1:length(TEsequence)) {
  seq_length[[f]] <- summary(TEsequence[[f]])$length
}  

seq_length <- as.data.frame(seq_length)  
seq_length$repeat_id <- attr(TEsequence, "name")  
seq_length$repeat_id <- gsub('\\#.*','', seq_length$repeat_id) #tidies up repbase naming conventions  

#merge   
test <- merge(noisonosim, seq_length, by = "repeat_id")   
test <- dplyr::rename(test, reference_length = seq_length)  
test <- test %>%  mutate(., Query_Length = ((qry_end - qry_start) + 1)) #overall length of fragment  
test <- test %>% group_by(repeat_id, matching_repeat, qry_id) %>% mutate(., lowextremety = min(qry_start)) #the earliest the fragment appears in the scaffold  
test <- test %>% group_by(repeat_id, qry_id, matching_repeat) %>% mutate(., highextremety = max(qry_end))  #the latest the fragment appears in the scaffold  


test <- mutate(test, mergedfraglength = ((highextremety - lowextremety)+1)) #overall merged fragment length  


test$length_check <- ifelse(test$mergedfraglength <= test$reference_length, "YES", "NO")  

test1 <- filter(test, length_check == "YES") #put all cases where the merged fragment length is less than the refbase entry into one database  
test2 <- filter(test, length_check == "NO") #put all cases where the merged fragment length is longer than the refbase entry into one database  

#####test2#####  
#next snippet of code keeps merged elements that should have been merged seperate  



test2 <- test2 %>% dplyr::select(-c('lowextremety', 'highextremety', 'mergedfraglength'))  

test2 <- rename(test2, lowextremety = qry_start)  
test2 <- rename(test2, highextremety = qry_end)  
test2 <- rename(test2, mergedfraglength = qry_width)  

##########################  


noisonosim <- rbind(test1, test2) #merge  


noisonosim <- noisonosim %>% dplyr::select(-c(Query_Length, qry_width, in_repeat_start, in_repeat_end, no_bp_in_complement, qry_start, qry_end, qry_left, sw_score)) #remove all these columns  
noisonosim$perc_del <- as.numeric(noisonosim$perc_del) #ensure these column are numeric  
noisonosim$perc_div <- as.numeric(noisonosim$perc_div)  
noisonosim$perc_insert <- as.numeric(noisonosim$perc_insert)  

#calculate a new mean perc del across the novel merged fragments copy   
noisonosim <- noisonosim %>% group_by(qry_id, repeat_id, matching_repeat, mergedfraglength) %>% mutate(., Copy_perc_div = mean(perc_div))   
noisonosim <- noisonosim %>% group_by(qry_id, repeat_id, matching_repeat, mergedfraglength) %>% mutate(., Copy_perc_del = mean(perc_del))  
noisonosim <- noisonosim %>% group_by(qry_id, repeat_id, matching_repeat, mergedfraglength) %>% mutate(., Copy_perc_insert = mean(perc_insert))  


#remove old 'perc del' values  
noisonosim <- noisonosim %>% dplyr::select(-c(perc_div, perc_del, perc_insert))  

noisonosim <- dplyr::rename(noisonosim, perc_div = Copy_perc_div) #rename  
noisonosim <- dplyr::rename(noisonosim, perc_del = Copy_perc_del)  
noisonosim <- dplyr::rename(noisonosim, perc_insert = Copy_perc_insert)  


noisonosim <- unique(noisonosim) #this step ensures that what was two fragments now merged is only represented once  



#3)Remove TEs found in different isoforms of the same "gene"   
#########NOTE - this step is not needed when working with genome data##############  


noisonosim$Gene <- noisonosim$qry_id #make duplicated transcript file  

noisonosim$Gene <- sub("_[^_]+$","", noisonosim$Gene) #removes isoform identifier  

noisonosim$isoform <- noisonosim$qry_id #make another duplicated transcript file  

noisonosim$isoform <- gsub(".*i","", noisonosim$isoform) #removes isoform identifier  


#for a given TE how many isoforms is it found on  
noisonosim <- noisonosim %>% group_by(Gene, repeat_id, matching_repeat) %>% mutate(isoform_number = n_distinct(isoform))   

check <- filter(noisonosim, isoform_number == 1) #if TE is only found in one isoform then no worry  



#What about TEs found in multiple isoforms?  
check2 <- filter(noisonosim, isoform_number > 1)  
#only keep isoform which has the most TE content  
check2 <- check2 %>% group_by(repeat_id, qry_id, matching_repeat) %>% mutate(totalTElength_pertranscript = sum(mergedfraglength)) #total TE content per element per transcript (ie if fragments were merged)  
check2 <- check2 %>% ungroup() %>% group_by(repeat_id, Gene, matching_repeat) %>% filter(.,totalTElength_pertranscript ==max(totalTElength_pertranscript))  

#however - this will correctly keep two TE fragments on the same isoform but incorrectly keep TEs on different isoforms if they have the same merged length (common)  
check2 <- check2 %>% group_by(Gene, repeat_id, matching_repeat) %>% mutate(isoform_number_two = n_distinct(isoform)) #are TEs on same isoform or not?  

check3 <- check2 %>% ungroup() %>% filter(isoform_number_two > 1) #check 3 contains those incorrect - found on two different isoforms  
check2 <- check2 %>% filter(isoform_number_two == 1) #check 2 contains those correct - fragments found on same isoform  
check2 <- check2 %>% dplyr::select(-c(isoform_number_two, totalTElength_pertranscript)) #tidy  




check3 <- check3 %>% distinct(Gene, repeat_id, matching_repeat, .keep_all = TRUE) #only keep one instance (at random) if TEs on different isoforms have same merged frag length  
#distinct (only keep one entry) if sample, gene, repeat and matching repeat all identical  

check3 <- check3 %>% dplyr::select(-c(isoform_number_two, totalTElength_pertranscript)) #tidy  


#time to merge those TEs found in only one isoform (check), largest TEs containing isoform if found in multiple isoforms  (keeping both fragments on same isoform if necessary)  (check2), or if no largest TE containing isoform one entry kept at random  
check <-as.data.frame(check)  
check2 <- as.data.frame(check2)  
check3 <- as.data.frame(check3)  
finalcheck <- rbind(check, check2, check3)  

noisonosim <- as.data.frame(finalcheck)  


#4) remove fragments less than 80bp  
noisonosim <- filter(noisonosim, mergedfraglength >= 80)  

#5) Write output file  
#quick clean up  
noisonosim <- noisonosim %>% dplyr::select(-c(length_check, isoform_number))  
noisonosim <- rename(noisonosim, merged_qrystart = lowextremety)  
noisonosim <- rename(noisonosim, merged_qryend = highextremety)  


setwd(i)  
write.csv(noisonosim, file = paste0(j, "_RM_TRIPS.csv"))  


### RM_TRIPS_END


## Step7 adjust RM Trips output to divsum-file format

###with R using R version 4.1.2 (2021-11-01)  
###gtools_3.9.4; stringr_1.5.0; tidyverse_2.0.0; reshape2_1.4.4; ggplot2_3.4.1; dplyr_1.1.1  
###I made this myself, so please don't hate if something does not work


############################################
### RMTrips output to divsum-file format ###
########## Viktoria Bednarski ##############
############################################


rm(list = ls())  

library(reshape2)  
library(dplyr)  
library(tidyverse)  
library(gtools)  
library(stringr)  
library(ggplot2)  

rm(list = ls())  

RM_Trips_table <- read.csv('/home/jens/Viki/stuff_for_RM_Trips/Ppr germany/alt2/tidy_noasterisk_alt2_Chr9_1_1.fasta.out_RM_TRIPS.csv')  

#Reducing the RMTrips output to the required factors  
RMTrips_small <- RM_Trips_table[, c("matching_class", "mergedfraglength", "perc_div")]  

#Group by matching_class and perc_div, and then sum the values in mergedfraglength  
merged_table <- RMTrips_small %>%  
  group_by(matching_class, perc_div) %>%  
  summarise(sum_mergedfraglength = sum(mergedfraglength))  

######################this is optional, just to check if perc_div is numeric##################################  

#Check the data types of entries in merged_table  
data_types <- sapply(merged_table, class)  

#Print the data types  
print(data_types)  

##############################################################################################################  

#Define the intervals for perc_div            
##########this might need some refinement, because it doesn't group values of 0,00 into the (0,1] interval and i don't know why, so this is why they have to be manually renamed in a later step and then added  
merged_table <- merged_table %>%  
  mutate(interval = cut(perc_div, breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80)))  


#Group by matching_class, interval, and sum the values in sum_mergedfraglength  
interval_merged_table <- merged_table %>%  
  group_by(matching_class, interval) %>%  
  summarise(sum_sum_mergedfraglength = sum(sum_mergedfraglength))  

#Replace NA entries in interval with a suitable identifier  
interval_merged_table$interval[is.na(interval_merged_table$interval)] <- "(0,1]"  

#to change the names of the TE types  
names <- read.csv('/home/jens/Desktop/converted table/TE types names.csv')  

#Perform a left join between dataframe1 and dataframe2 based on the "names" column  
interval_merged_table_nn <- merge(interval_merged_table, names, by.x = "matching_class", by.y = "names", all.x = TRUE)  

#Replace the values in column1 with the corresponding values from the "new_names" column  
interval_merged_table_nn$matching_class <- interval_merged_table_nn$new_names  

#Remove the "new_names" column if it's no longer needed  
interval_merged_table_nn <- subset(interval_merged_table_nn, select = -c(new_names))  


#Pivot the table to convert matching_class values as column names and perc_div (now interval) values as row names  
converted_table <- interval_merged_table_nn %>%  
  pivot_wider(names_from = matching_class, values_from = sum_sum_mergedfraglength, values_fn = sum)  

#To avoid annoying error message, convert the tibble to a data frame  
converted_table <- as.data.frame(converted_table)  

#Sort the table by the interval column in ascending order  
sorted_data <- converted_table[order(converted_table$interval), ]  

#Replace interval values with ascending numbering  
sorted_data <- sorted_data %>%   
  mutate(interval = seq_along(interval) - 1)  

#Replace NA entries with 0  
sorted_data <- replace(sorted_data, is.na(sorted_data), 0)  

#Assign a new name to the interval column, so it works with the other script  
colnames(sorted_data)[1] <- "Div"  



#Specify the desired column order  
new_order <- c("Div","DTA_hAT","DTA_hAT_MITE","DTA_hAT_unknown","DTB_PiggyBac","DTB_PiggyBac_MITE","DTC_CACTA","DTC_CACTA_MITE","DTC_CACTA_unknown","DTH_Harbinger","DTH_Harbinger_MITE","DTM_Mutator","DTM_Mutator_MITE","DTM_Mutator_unknown","DTP_P","DTT_TC1_Mariner","DTT_TC1_Mariner_MITE","DTT_TC1_Mariner_unknown","DXX","DHH_Helitron","RLB_Bel_Pao","RLC_Copia","RLE_ERV","RLG_Gypsy","RLX","RII_I","RIJ_Jockey","RIL_L1","RIR_R2","RIT_RTE","RIX","RPP_Penelope","RST_tRNA","RXX_nLTR","RXX","XXX")  

#Retain only the common elements between new_order and existing column names  
cols_to_reorder <- intersect(new_order, colnames(sorted_data))  

#Reorder the columns of sorted_data  
sorted_data_reordered <- sorted_data[, cols_to_reorder]  

#Set the working directory  
setwd("/home/jens/Viki/stuff_for_RM_Trips/Ppr germany/alt2")  

#Save the converted table as CSV  
write.csv(sorted_data_reordered, file = "Ppr_germany.alt2.Chr9.fasta.csv", row.names = FALSE)  


## Step8 Plot repeat divergence landscapes
###with R using R version 4.1.2 (2021-11-01)  
###taken and adjusted from KristinaGagalova's comment: https://github.com/oushujun/EDTA/issues/92  
###watch out if you are using this script, then you should not go back to step6 in the same R window, because there is some kind of package conflict and then the script of step6 is not working properly and giving error messages which leads to undesirable and weird results in the repeat divergence landscapes, working in two seperate windows seems to work though  

#############################
### Plot Kimura distance ####
#############################

rm(list = ls())  


library(reshape)  
library(ggplot2)  
library(hrbrthemes)  
library(tidyverse)  
library(gridExtra)  

KimuraDistance <- read.csv("/home/jens/Viki/stuff_for_RM_Trips/Ppr germany/alt2/Ppr_germany.alt2.Chr9.fasta.csv",sep=",")  

#add here the genome size in bp, remember to change accordingly  
genomes_size=16549588  

kd_melt = melt(KimuraDistance,id="Div")  
kd_melt$norm = kd_melt$value/genomes_size * 100  

#Define a named vector of custom colors for specific items  
custom_colors <- c("DHH_Helitron" = "#91A6FF", "DTA_hAT" = "#F76CC6", "DTA_hAT_MITE" = "#F76CC6", "DTA_hAT_unknown" = "#F76CC6", "DTB_PiggyBac" = "#E85DBA", "DTB_PiggyBac_MITE" = "#E85DBA", "DTC_CACTA" = "#DA4DAD", "DTC_CACTA_MITE" = "#DA4DAD", "DTC_CACTA_unknown" = "#DA4DAD", "DTH_Harbinger" = "#CB3EA1", "DTH_Harbinger_MITE" = "#CB3EA1", "DTM_Mutator" = "#BC2E94", "DTM_Mutator_MITE" = "#BC2E94", "DTM_Mutator_unknown" = "#BC2E94", "DTP_P" = "#AD1F88", "DTT_TC1_Mariner" = "#9F0F7B", "DTT_TC1_Mariner_MITE" = "#9F0F7B", "DTT_TC1_Mariner_unknown" = "#9F0F7B", "DXX" = "#90006F", "RLB_Bel_Pao" = "#1EFFBC", "RLC_Copia" = "#17DEA3", "RLE_ERV" = "#0FBD89", "RLG_Gypsy" = "#089B70", "RLX" = "#007A56", "RII_I" = "#5EAE5F", "RIJ_Jockey" = "#4D9F4E", "RIL_L1" = "#448D45", "RIR_R2" = "#3A7B3B", "RIT_RTE" = "#2F6930", "RIX" = "#235724", "RPP_Penelope" = "#32746D", "RST_tRNA" = "#2A6B67", "RXX_nLTR" = "#216261", "RXX" = "#104F55", "XXX" = "grey")  


#Plot with custom colors and legend  
ggplot(kd_melt, aes(fill = variable, y = norm, x = Div)) +
  geom_bar(position = "stack", stat = "identity", color="black") +
  scale_fill_manual(values = custom_colors, na.value = "blue") +
  theme_classic() +
  xlab("Kimura substitution level") +
  ylab("Percent of the genome") +
  coord_cartesian(xlim = c(0, 40), ylim = c(0, 2)) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 12)) +
  labs(fill = "TE Types")  

