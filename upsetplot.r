#topGo script 
library("UpSetR")

#Set path
setwd("C:/Users/vicme/OneDrive/Desktop/PhD/pawpaw_project/Finalizing pipeline")

orthogroups_df<- read.table("Orthogroups.GeneCount.tsv",  header=T, stringsAsFactors = F)

#All species
selected_species <- colnames(orthogroups_df)[2:(ncol(orthogroups_df) -1)] 
selected_species
ncol(orthogroups_df)

#Need to take out list of orthogroups that didnt have any results from interproscan, retroelements, or transposon/transposase related
include_ogs <- read.table("include_orthogroups.txt", header = FALSE, stringsAsFactors = FALSE)[,1]

# Filter dataframe
filtered_df <- orthogroups_df[orthogroups_df$Orthogroup %in% include_ogs, ]


filtered_df[filtered_df > 0] <- 1
upset(filtered_df, nsets = ncol(filtered_df), sets = rev(selected_species), 
      keep.order = T, order.by = "freq", number.angles = 30, empty.intersections = "on")
