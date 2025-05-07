rm(list=ls())
library(tidyverse)

#(Pseudopanurgus,Peponapis,Tetraloniella,Syntrichalonia,Cemolobus,Micralictoides)
#Eucerini treated according to https://academic.oup.com/isd/article/7/4/3/7222693?login=true 
#Pseudopanurgus https://resjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/syen.12530
#Micralictoides treated as subgenus of Duforea
#Dufourea subg. Micralictoides Timberlake, 1939. Ent. Soc. Amer., Ann. 32: 397.

# Read the data
all_genera <- read.table('../data/all-genera.txt', header = TRUE)

# Extract unique rows and rename the column to 'genus'
unique_genus <- unique(all_genera)
colnames(unique_genus) <- "genus"

# Write to CSV
write.csv(unique_genus, "../data/unique_genus.csv", row.names = FALSE, , quote = FALSE)


# Perform multiple generic replacements<- did not use this part of the code
replace_species <- function(x) {
  x <- gsub("Pseudopanurgus", "Protandrena", x, ignore.case = TRUE)
  x <- gsub("Peponapis", "Xenoglossa", x, ignore.case = TRUE)
  x <- gsub("Tetraloniella", "Xenoglossa", x, ignore.case = TRUE)
  x <- gsub("Syntrichalonia", "Xenoglossa", x, ignore.case = TRUE)
  x <- gsub("Cemolobus", "Xenoglossa", x, ignore.case = TRUE)
  x <- gsub("Micralictoides", "Dufourea", x, ignore.case = TRUE)
  return(x)
}

# Apply the replacement function to all columns in the dataframe for updated nomenclature
unique_genera <- unique_genera %>%
  mutate_all(~ replace_species(.))

#write updated dataframe
write_csv(unique_genera,'data/genera.csv')

