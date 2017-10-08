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
if (max(width(seqs)) < 10000) {
  aln <- msa(seqs, method="Muscle")
  #print(aln, show="complete")
  writeXStringSet(as(aln, "DNAStringSet"), 
                  filepath=paste0(name, ".aln"))
  con <- consensusString(consensusMatrix(aln), threshold=1)
  write(con, file=paste0(name, ".consensus"))
} else {
  write("Error: sequences too long. Alignment in R using Muscle gives fatal error...",
        file=paste0(name, ".consensus"))
}



