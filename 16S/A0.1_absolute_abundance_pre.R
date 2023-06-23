#!/usr/bin/env Rscript

arg <- commandArgs(T)
if(length(arg) > 4){
        cat("Argument: Out_Dir Data_File\n")
        quit('no')
}
print(arg)

setwd(arg[1])

merge.file <- arg[2]    # "merged.dmp"
asvclass.file <- arg[3] # "Total.asvclass.txt"
dbpy.file <- arg[4]     # "dbpy.txt"
library(tidyverse)

merge.dta <- read.delim(
  merge.file, sep = "|",
  col.names = c("spid", "new", "dd")
) %>% .[-3]

data <- read_tsv(asvclass.file) %>% 
  left_join(., merge.dta)

data$new[is.na(data$new)] <- data$spid[is.na(data$new)]

data <- data[-2]
colnames(data)[2] <- "spid"

taxonomy.dta <- read_tsv(dbpy.file)


Data <- left_join(data, taxonomy.dta) %>% 
  filter(!is.na(taxonomy)) %>% 
  group_by(sample, taxonomy) %>% 
  summarise(counts = n()) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = sample,
    values_from = counts,
    values_fill = 0
  ) %>% 
  separate(
    taxonomy,
    sep = "\\|",
    into = c("kindom", "phylum", "class", "order",
             "family", "genus", "species")
  )

neibiao <- colSums(Data[str_detect(Data$family, pattern = "f__Planococcaceae"),8:ncol(Data)])
Data[8:ncol(Data)] <- Data[8:ncol(Data)]/(unname(neibiao)*100)

write_tsv(Data, "Total.taxonomy.absolute.abundance.txt")

Data <- Data[-which(str_detect(Data$phylum, pattern = "p__unclassified") | str_detect(Data$genus, pattern = "g__unclassified")),-c(1, 3:5)]
Data$phylum <- str_remove(Data$phylum, pattern = "p__")
Data$phylum <- str_replace_all(Data$phylum, " ", "_")
Data$genus <- str_remove(Data$genus, pattern = "g__")
Data$genus <- str_replace_all(Data$genus, " ", "_")
Data$species <- str_remove(Data$species, pattern = "s__")
Data$species <- str_replace_all(Data$species, " ", "_")

Phylum.dta <- Data[-c(2, 3)] %>% 
  rename(ID = phylum) %>% 
  write_tsv("Total.Phylum.absolute.abundance.txt")

Genus.dta <- Data[-c(1, 3)] %>% 
  rename(ID = genus) %>% 
  write_tsv("Total.Genus.absolute.abundance.txt")

Species.dta <- Data[-c(1, 2)] %>% 
  rename(ID = species) %>% 
  write_tsv("Total.Species.absolute.abundance.txt")

Plygen.dta <- Data %>% 
  pivot_longer(
    colnames(Data)[4:ncol(Data)],
    names_to = "sample",
    values_to = "abundance"
  ) %>% 
  group_by(phylum, genus, sample) %>% 
  summarise(abundance = sum(abundance)) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = "sample",
    values_from = "abundance"
  )

Phylum.dta <- Plygen.dta[-c(2)] %>% 
  rename(ID = phylum) %>% 
  write_tsv("Total.Phylum.absolute.abundance.txt")
Genus.dta <- Plygen.dta[-c(1)] %>% 
  rename(ID = genus) %>% 
  write_tsv("Total.Genus.absolute.abundance.txt")

Plygen.dta$phylum <- paste(Plygen.dta$phylum, Plygen.dta$genus, sep = " ")
Plygen.dta2 <- Plygen.dta[-c(2)] %>% 
  rename("Phyla Genera" = phylum) %>% 
  write_tsv("Total.PlyGen.absolute.abundance.txt")

Plyspe.dta <- Data
Plyspe.dta$phylum <- paste(Plyspe.dta$phylum, Plyspe.dta$species, sep = " ")
Plyspe.dta2 <- Plyspe.dta[-c(2, 3)] %>% 
  rename("Phyla Species" = phylum) %>% 
  write_tsv("Total.PlySpe.absolute.abundance.txt")
