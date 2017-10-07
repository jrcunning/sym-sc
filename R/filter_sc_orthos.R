library(stringr)
library(tidyverse)

#args = commandArgs(trailingOnly=TRUE)

dat <- read.delim("fastOrtho_out/ortho.table",header=F,sep="\t")
colnames(dat) <- c("Group", "Genes", "Taxa", "IDs")
dat.all <- filter(dat, Taxa>=4)
dim(dat.all) #7430 ortholog groups occur in all 3 sym genomes

counts <- data.frame(matrix(nrow=nrow(dat.all),ncol=4))
rownames(counts) <- dat.all$Group
colnames(counts) <- c("Smic","SymB","Sgor","Skaw")

counts$Smic <- apply(dat.all,1,function(x) str_count(x[4],"Smic.genome.annotation.pep.longest"))
counts$SymB <- apply(dat.all,1,function(x) str_count(x[4],"symbB.v1.2.augustus.prot"))
counts$Sgor <- apply(dat.all,1,function(x) str_count(x[4],"SymbC1.Gene_Models.PEP.fasta"))
counts$Skaw <- apply(dat.all,1,function(x) str_count(x[4],"SymbF.Gene_Models.PEP.fasta"))

# Find groups with a single occurrence in each genome (purported single copy loci)
sc <- counts[which(rowSums(counts)==4), ]
nrow(sc)
sc.fams <- rownames(sc)

write.table(sc.fams, file="sc_orthos/sc.names.txt", row.names=F, col.names=F, quote=F)

