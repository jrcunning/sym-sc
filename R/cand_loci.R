library(Biostrings)
library(stringr)

for (consensus in list.files(path="sc_orthos/each", pattern="*.consensus", full.names=T)) {
  name <- basename(tools::file_path_sans_ext(consensus))
  # Read in consensus sequence
  cons <- readLines(consensus)
  # Replace question marks with dashes
  cons <- str_replace_all(cons, "\\?", "-")
  # Convert consensus sequence to DNAStringSet object
  cons <- tryCatch({
    as(cons, "DNAStringSet")
  }, error = function(e) {
    return(NULL)
  })
  if (!is.null(cons)) {
    # Look for instances of 25 consecutive bases (candidate primer regions)
    pat <- DNAString(paste0(rep("N", 20), collapse=""))
    cand <- vmatchPattern(pat, cons, fixed=F, max.mismatch=0)
    # If there are more than one conserved runs > 25
    if (length(cand[[1]]) >=2) {
      # Check that there is a potential amplicon >40 and <150 bp (which includes length of one primer)
      amp.lengths <- apply(combn(start(cand)[[1]], 2), 2, diff)
      good.cand <- any(amp.lengths > 50 & amp.lengths < 2000)
      # If there are good candidate primer sites, write to file
      if (good.cand==T) {
        region <- subseq(cons, min(start(cand)), max(end(cand)))
        names(region) <- name
        writeXStringSet(region, filepath=paste0("cand_loci/", name, ".region"))
      }
    }
  }
}



  
