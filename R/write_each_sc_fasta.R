library(Biostrings)

# Read in all single copy sequences from all Symbiodinium
allsym <- readDNAStringSet("sc_orthos/allsym_sc.fasta")
names(allsym) <- lapply(strsplit(names(allsym), split=" "), "[[", 1)

# Read in groupings for each single copy group
orthos <- read.table("sc_orthos/sc.ortho.table.trimmed", stringsAsFactors = F)

# Separate sequences into individual groups and write to fasta
for (i in 1:nrow(orthos)) {
  group <- orthos[i, 1]
  names <- orthos[i, 2:5]
  seqs  <- allsym[names(allsym) %in% names]
  writeXStringSet(seqs, filepath=paste0("sc_orthos/each/", group, ".fasta"))
}
