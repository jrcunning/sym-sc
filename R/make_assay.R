library(Biostrings)
library(stringr)

loci <- readDNAStringSet("cand_loci/all_candidates.fasta")

runs <- str_extract_all(loci, "[^-]*[^-]")
primers <- lapply(runs, function(x) x[which(nchar(x)>=20)])

assays <- data.frame(
  assay=names(loci),
  forward=unlist(lapply(primers, "[[", 1)),
  reverse=as.character(reverseComplement(as(unlist(lapply(primers, "[[", 2)), "DNAStringSet")))
)

write.table(assays, file="cand_loci/pcr/assays.txt", row.names=F, col.names=F, quote=F, sep=" ")
