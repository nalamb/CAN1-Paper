
##### Load required packages #####
require ("lattice")
library(ggplot2)
library(gplots)
library(latticeExtra)
library(dplyr)
library(pheatmap)
library(gridExtra)
library(grid)
library(reshape)
library(devtools) 
library(tidyr)
library(corrplot)
library(Hmisc)
library(xlsx)
library(RColorBrewer)

##### Permissive filter #####
dt <- read.csv("AllVariantsEdit.txt",sep = "\t")
# Set up a list of all permissive groupings
permissiveGroups = c("WT_perm",
                     "rnr1Y285F_perm",
                     "WT_2424_perm",
                     "rnr1Y285A_perm",
                     "rnr1D57N_perm")
# Calculate out the total background frequency of permissive samples
permissiveAvgFreq <- mean(dt[dt$Group %in% permissiveGroups,]$Frequency)
AvgFreq <- mean(dt$Frequency)
#create a dataframe with columns to store the filtered data in
filter <- data.frame(list(Reference.Position = 0, numSamps=0,avgFreq=0,bpFilt = 0))
# Rules for filtering
# for each basepair in the CAN1 range -- iterate through and look for various things :)
for (i in 31639:33479){
  #calculate the total number of variants in the permissive samples for that basepair
  numSamps = nrow(dt[(dt$Reference.Position %in% c(i) & dt$Group %in%permissiveGroups),])
  #calculate the maximum frequency for a given variant at a specific basepair
  maxFreq = max(dt[(dt$Reference.Position %in% c(i) & dt$Group %in%permissiveGroups),]$Frequency)
  #calculate the average frequency in the permissive samples for a given basepair
  avgFreq = mean(dt[(dt$Reference.Position %in% c(i) & dt$Group %in%permissiveGroups),]$Frequency)
  
  #if there are no variants present at a given basepair
  if (numSamps < 1){
    #set the cutoff to the total background frequency for all permissive samples
    bpFilt = permissiveAvgFreq
    # else if there are variants present
  } else if (numSamps > 0){
    #set the cutoff to the maximum frequency of a permissive variant
    bpFilt = avgFreq
  }
  
  # Lastly, if the new cutoff is less than the average background of all perm variants
  if (bpFilt < permissiveAvgFreq){
    #set it to the avg perm variant frequency (.108)
    bpFilt = permissiveAvgFreq
  }
  
  # finally -- create a new row dataframe with this information for a given basepair
  newRow <- data.frame(list(Reference.Position = i,
                            numSamps = numSamps,
                            avgFreq = avgFreq,
                            bpFilt = bpFilt
  ))
  ## and append it on to our filter data.frame
  filter <- rbind(filter,newRow)
}
# Once the filter dataframe has been created, apply the filter to our original dataframe
dt.filt <- left_join(dt,filter,by="Reference.Position") %>% filter(Frequency > bpFilt)
#Subset on quality, growth condition and to include only 1 technical replicate
#Galactose induced samples were not included in this manuscript
dt.filt <- dt.filt %>% filter(Average.quality >30)
Allvariants_table_notech <- dt.filt %>% filter(Tech_rep == "FALSE")
Allvariants <- Allvariants_table_notech %>% filter(GAL_status == "FALSE")
Allvariants <- Allvariants %>% filter(Growth.Condition == "CAN")


##### Mutation spectra are reported as normalized frequency (counts) or total frequency #####

##To plot percentages of insertions, deletions and different SNVs, plotting total variant frequency
All_clustID_grouped_sum <- aggregate(Allvariants$Frequency, 
                                     by =list(Allvariants$Group, Allvariants$Clust.I.D.), FUN=sum)
##To plot percentages of insertions, deletions and different SNVs, plotting counts of unique variants
#All_clustID_grouped_length <- aggregate(Allvariants$Frequency, 
#                                     by =list(Allvariants$Group, Allvariants$Clust.I.D.), FUN=length)

#Shifting the table structure 
Allvariants_table <- All_clustID_grouped_sum %>% select(Group.2)
Allvariants_unique <- unique(Allvariants_table)
uniqueList = (unique(Allvariants$Group))
uniqueTable <- unique(Allvariants$Group)
for (n in uniqueList){
  tempTable <- All_clustID_grouped_sum %>% filter(Group.1 == n) 
  tempTableNew <- tempTable %>% select(Group.2, x)
  names(tempTableNew)[2]<-n
  Allvariants_unique <- full_join(Allvariants_unique, tempTableNew, by = "Group.2")
}
Allvariants_unique <- Allvariants_unique %>% replace(is.na(.),0.0)
Allvariants_unique_stripped = subset(Allvariants_unique, select = -c(Group.2))
NewTableToMakePlots <- subset(Allvariants_unique, select = c(Group.2))
Allvariants_unique.matrix <- as.matrix(Allvariants_unique_stripped)
NEW_Allvariants_unique <- sweep(Allvariants_unique.matrix,MARGIN=2,FUN="/",STATS=colSums(Allvariants_unique.matrix))
NewDF <- as.data.frame(NEW_Allvariants_unique)
NewTableToMakePlots2 <- cbind(NewTableToMakePlots, NewDF)

#converting to long format to be able to make barcharts
data_long <- gather(NewTableToMakePlots2, genotype, frequency, msh2:WT, factor_key=TRUE)
#SNV, deletion and insertion subsets
data_long_SNVs <- data_long %>% filter(data_long$Group.2 %in% c("CG>TA", "CG>AT", "CG>GC", "TA>AT", "TA>CG", "TA>GC" ))
data_long_deletions <- data_long %>% filter(data_long$Group.2 %in% c("A/T-1", "G/C-1", ">1bpdel"))
data_long_insertions <- data_long %>% filter(data_long$Group.2 %in% c("A/T+1", "G/C+1", ">1bp-ins"))

#Example barplots
ggplot(data_long_SNVs[data_long_SNVs$genotype %in% c("WT", "msh2", "msh3", "msh6", "rnr1D57N", "rnr1D57N_msh2", "rnr1D57N_msh3", "rnr1D57N_msh6", "rnr1Y285F", "rnr1Y285F_msh2", "rnr1Y285F_msh3", "rnr1Y285F_msh6", "rnr1Y285A", "rnr1Y285A_msh2", "rnr1Y285A_msh3", "rnr1Y285A_msh6"),], 
       aes(fill=Group.2, y=frequency, x=genotype)) + 
  geom_bar( colour = "black", stat="identity",position="stack") +
  #scale_fill_brewer(palette="Paired") +
  #scale_fill_manual(values=c("seashell3", "midnightblue", "plum3", "palevioletred4", "lightseagreen")) 
  ylab("Normalized Frequency") + theme(axis.title.x = element_blank(), legend.title = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5) ) + 
  scale_x_discrete(limits=c("WT", "msh2", "msh3", "msh6", "rnr1D57N", "rnr1D57N_msh2", "rnr1D57N_msh3", "rnr1D57N_msh6", "rnr1Y285F", "rnr1Y285F_msh2", "rnr1Y285F_msh3", "rnr1Y285F_msh6", "rnr1Y285A", "rnr1Y285A_msh2", "rnr1Y285A_msh3", "rnr1Y285A_msh6")) + 
  theme(axis.text.x = element_text(size = 6)) 

ggplot(data_long_deletions[data_long_deletions$genotype %in% c("WT", "msh2", "msh3", "msh6", "rnr1D57N", "rnr1D57N_msh2", "rnr1D57N_msh3", "rnr1D57N_msh6", "rnr1Y285F", "rnr1Y285F_msh2", "rnr1Y285F_msh3", "rnr1Y285F_msh6", "rnr1Y285A", "rnr1Y285A_msh2", "rnr1Y285A_msh3", "rnr1Y285A_msh6"),],
       aes(fill=Group.2, y=frequency, x=genotype)) + 
  geom_bar( colour = "black", stat="identity",position="stack") +
  #scale_fill_brewer(palette="Paired") +
  #scale_fill_manual(values=c("seashell3", "midnightblue", "plum3", "palevioletred4", "lightseagreen")) 
  ylab("Normalized Frequency") + theme(axis.title.x = element_blank(), legend.title = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  scale_x_discrete(limits=c("WT", "msh2", "msh3", "msh6", "rnr1D57N", "rnr1D57N_msh2", "rnr1D57N_msh3", "rnr1D57N_msh6", "rnr1Y285F", "rnr1Y285F_msh2", "rnr1Y285F_msh3", "rnr1Y285F_msh6", "rnr1Y285A", "rnr1Y285A_msh2", "rnr1Y285A_msh3", "rnr1Y285A_msh6")) + theme(axis.text.x = element_text(size = 6)) 

ggplot(data_long_insertions[data_long_insertions$genotype %in% c("WT", "msh2", "msh3", "msh6", "rnr1D57N", "rnr1D57N_msh2", "rnr1D57N_msh3", "rnr1D57N_msh6", "rnr1Y285F", "rnr1Y285F_msh2", "rnr1Y285F_msh3", "rnr1Y285F_msh6", "rnr1Y285A", "rnr1Y285A_msh2", "rnr1Y285A_msh3", "rnr1Y285A_msh6"),],
       #ggplot(data_long_insertions[data_long_insertions$genotype %in% c("WT", "msh2", "msh3", "msh6", "rnr1D57N", "rnr1D57N_msh2", "rnr1D57N_msh3", "rnr1D57N_msh6"),],
       aes(fill=Group.2, y=frequency, x=genotype)) + 
  geom_bar( colour = "black", stat="identity",position="stack") + 
  #scale_fill_brewer(palette="Paired") +
  #scale_fill_manual(values=c("seashell3", "midnightblue", "plum3", "palevioletred4", "lightseagreen")) 
  ylab("Normalized Frequency") + theme(axis.title.x = element_blank(), legend.title = element_blank(), axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  scale_x_discrete(limits=c("WT", "msh2", "msh3", "msh6", "rnr1D57N", "rnr1D57N_msh2", "rnr1D57N_msh3", "rnr1D57N_msh6", "rnr1Y285F", "rnr1Y285F_msh2", "rnr1Y285F_msh3", "rnr1Y285F_msh6", "rnr1Y285A", "rnr1Y285A_msh2", "rnr1Y285A_msh3", "rnr1Y285A_msh6")) + theme(axis.text.x = element_text(size = 6),) 


##### COSMIC Analysis #####
#COSMIC_sbs_GRC38_table <- read.table("COSMIC_v3.2_SBS_GRCh38.txt", header=TRUE)
##Read in the normalized frequencies for our dataset combined our SNV tri context data
dt.COSMIC <- read.csv("COSMIC_CAN1.csv")

##Heatmap and get correlation values
dt.COSMIC_stripped = subset(dt.COSMIC, select = -c(SNP.ID))
dt.COSMIC.matrix <- as.matrix(dt.COSMIC_stripped)
dt.COSMIC.cor = cor(dt.COSMIC.matrix, method= c("spearman"))
corrplot(dt.COSMIC.cor, type="upper")
corrplot(dt.COSMIC.cor)

palette = colorRampPalette(c("violetred", "mediumpurple", "lightblue")) (60)
heatmap.2(x = dt.COSMIC.cor, col = palette, symm = TRUE, key=TRUE, trace="none", margins = c(10,10))
write.table(dt.COSMIC.cor, file="SpearmanCorrelationSUM.txt", sep="\t")


