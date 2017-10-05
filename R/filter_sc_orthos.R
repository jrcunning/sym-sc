library(stringr)
library(tidyverse)

#args = commandArgs(trailingOnly=TRUE)

dat <- read.delim("protein.ortho.table",header=F,sep="\t")
colnames(dat) <- c("Group", "Genes", "Taxa", "IDs")
dat.all <- filter(dat, Taxa>=3)
dim(dat.all) #7430 ortholog groups occur in all 3 sym genomes

counts <- data.frame(matrix(nrow=nrow(dat.all),ncol=3))
rownames(counts) <- dat.all$Group
colnames(counts) <- c("Smic","SymB","Skaw")

counts$Smic <- apply(dat.all,1,function(x) str_count(x[4],"Smic.genome.annotation.pep.longest"))
counts$SymB <- apply(dat.all,1,function(x) str_count(x[4],"symbB.v1.2.augustus.prot"))
counts$Skaw <- apply(dat.all,1,function(x) str_count(x[4],"Symbiodinium_kawagutii.0819.final.gene"))

# Find groups with a single occurrence in each genome (purported single copy loci)
sc <- counts[which(rowSums(counts)==3), ]
nrow(sc)
sc.fams <- rownames(sc)

write.table(sc.fams, file="protein.sc.names.txt", row.names=F, col.names=F, quote=F)

