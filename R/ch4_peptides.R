RNA_nuc <- c("A", "C", "G", "U")

# generate codons
codons <- expand.grid(rep(list(RNA_nuc), 3), stringsAsFactors = FALSE)
codons <- do.call(paste0, as.list(codons[, 3:1]))

# amino acids
genetic_code <- c("K", "N", "K", "N", "T", "T", "T", "T", "R", "S", "R", "S", "I", "I", "M", "I", "Q", "H", "Q", "H", "P", "P", "P", "P", "R", "R", "R", "R", "L", "L", "L", "L", "E", "D", "E", "D", "A", "A", "A", "A", "G", "G", "G", "G", "V", "V", "V", "V", "*", "Y", "*", "Y", "S", "S", "S", "S", "*", "C", "W", "C", "L", "F", "L", "F") 
names(genetic_code) <- codons
TranslateCodon <- function(codon){
	translated_codon <- genetic_code[codon]
	names(translated_codon) <- NULL
	return(translated_codon)
}

TranslateRNA <- function(RNA){
	codons <- substring(RNA, seq(1, nchar(RNA), 3), seq(3, nchar(RNA), 3))
	peptide <- paste(TranslateCodon(codons), sep = "", collapse = "")
	n <- nchar(peptide)
	if (identical("*", substring(peptide, n))) {
		peptide <- substring(peptide, 1, n - 1)
	}
	return(peptide)
}
