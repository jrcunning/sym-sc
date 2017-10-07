library(Biostrings)
library(msa)

#for (fasta in list.files("sc_orthos/each", pattern="[0-9].fasta", full.names=T)) {
#  seqs <- readDNAStringSet(fasta)
#  name <- tools::file_path_sans_ext(fasta)
#  aln <- msa(seqs, method="Muscle")
#  #print(aln, show="complete")
#  writeXStringSet(as(aln, "DNAStringSet"), 
#                  filepath=paste0(name, ".aln"))
#  con <- consensusString(consensusMatrix(aln), threshold=1)
#  write(con, file=paste0(name, ".consensus"))
#}



args = commandArgs(trailingOnly=TRUE)

fasta <- args[1]
seqs <- readDNAStringSet(fasta)
name <- tools::file_path_sans_ext(fasta)
aln <- msa(seqs, method="Muscle")
#print(aln, show="complete")
writeXStringSet(as(aln, "DNAStringSet"), 
                filepath=paste0(name, ".aln"))
con <- consensusString(consensusMatrix(aln), threshold=1)
write(con, file=paste0(name, ".consensus"))


