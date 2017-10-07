all: sc_orthos/allsym_sc.fasta

sc_orthos/allsym_sc.fasta: sc_orthos/Skaw_sc.fasta sc_orthos/SymB_sc.fasta sc_orthos/Smic_sc.fasta sc_orthos/Sgor_sc.fasta
	cat $^ > $@

sc_orthos/Skaw_sc.fasta: sc_orthos/Skaw.sc.prots
	filter_fasta.py -f data/gene_models/SymbF.Gene_Models.CDS.fasta -o sc_orthos/Skaw_sc.fasta -s sc_orthos/Skaw.sc.prots

sc_orthos/SymB_sc.fasta: sc_orthos/SymB.sc.prots
	filter_fasta.py -f data/gene_models/symbB.v1.2.augustus.mrna.fa -o sc_orthos/SymB_sc.fasta -s sc_orthos/SymB.sc.prots

sc_orthos/Smic_sc.fasta: sc_orthos/Smic.sc.prots
	filter_fasta.py -f data/gene_models/Smic.genome.annotation.CDS.longest.sorted.fa -o sc_orthos/Smic_sc.fasta -s sc_orthos/Smic.sc.prots	

sc_orthos/Sgor_sc.fasta: sc_orthos/Sgor.sc.prots
	filter_fasta.py -f data/gene_models/SymbC1.Gene_Models.CDS.fasta -o sc_orthos/Sgor_sc.fasta -s sc_orthos/Sgor.sc.prots

sc_orthos/SymB.sc.prots: sc_orthos/sc.ortho.table.trimmed
	cut -f2 sc_orthos/sc.ortho.table.trimmed | cut -d' ' -f4 > sc_orthos/SymB.sc.prots

sc_orthos/Skaw.sc.prots: sc_orthos/sc.ortho.table.trimmed
	cut -f2 sc_orthos/sc.ortho.table.trimmed | cut -d' ' -f3 > sc_orthos/Skaw.sc.prots

sc_orthos/Sgor.sc.prots: sc_orthos/sc.ortho.table.trimmed
	cut -f2 sc_orthos/sc.ortho.table.trimmed | cut -d' ' -f2 > sc_orthos/Sgor.sc.prots

sc_orthos/Smic.sc.prots: sc_orthos/sc.ortho.table.trimmed
	cut -f2 sc_orthos/sc.ortho.table.trimmed | cut -d' ' -f1 > sc_orthos/Smic.sc.prots

sc_orthos/sc.ortho.table.trimmed: sc_orthos/sc.ortho.table
	cd sc_orthos && sed 's/(SymbF.Gene_Models.PEP)//; s/(SymbC1.Gene_Models.PEP)//; s/(Smic.genome.annotation.pep.longest)//; s/(symbB.v1.2.augustus.prot)//' sc.ortho.table | cut -f1,4 > sc.ortho.table.trimmed

sc_orthos/sc.ortho.table: sc_orthos/sc.names.txt
	grep -f sc_orthos/sc.names.txt fastOrtho_out/ortho.table > sc_orthos/sc.ortho.table

sc_orthos/sc.names.txt: R/filter_sc_orthos.R fastOrtho_out/ortho.table
	R --vanilla < R/filter_sc_orthos.R

fastOrtho_out/ortho.table:
	### Run FastOrtho on proteins
	FastOrtho \
	$$(printf "!!single_genome_fasta %s " data/proteins/* | tr "!" -) \
	--working_directory fastOrtho_out --project_name sym_ortho --blast_cpus 96

	##Format output file    
	sed -E "s/\s+\((\S+)\s+genes,(\S+)\s+taxa\)\:\s+/ \1 \2 /" fastOrtho_out/sym_ortho.end | sed -E -e 's/ /\t/' -e 's/ /\t/' -e 's/ /\t/' > fastOrtho_out/ortho.table
