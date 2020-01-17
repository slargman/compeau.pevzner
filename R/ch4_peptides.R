RNA_nuc <- c("A", "C", "G", "U")

# generate codons
codons <- expand.grid(rep(list(RNA_nuc), 3), stringsAsFactors = FALSE)
codons <- do.call(paste0, as.list(codons[, 3:1]))

# amino acids
genetic_code <- c("K", "N", "K", "N", "T", "T", "T", "T", "R", "S", "R", "S", "I", "I", "M", "I", "Q", "H", "Q", "H", "P", "P", "P", "P", "R", "R", "R", "R", "L", "L", "L", "L", "E", "D", "E", "D", "A", "A", "A", "A", "G", "G", "G", "G", "V", "V", "V", "V", "*", "Y", "*", "Y", "S", "S", "S", "S", "*", "C", "W", "C", "L", "F", "L", "F") 
names(genetic_code) <- codons

# amino acid masses
amino_acid_mass <- c(57, 71, 87, 97, 99, 101, 103, 113, 113, 114, 115, 128, 128, 129, 131, 137, 147, 156, 163, 186)
names(amino_acid_mass) <- c("G", "A", "S", "P", "V", "T", "C", "I", "L", "N", "D", "K", "Q", "E", "M", "H", "F", "R", "Y", "W")

TranslateCodon <- function(codon){
	translated_codon <- genetic_code[codon]
	names(translated_codon) <- NULL
	return(translated_codon)
}

TranslateRNA <- function(RNA, vector = FALSE){
	codons <- substring(RNA, seq(1, nchar(RNA), 3), seq(3, nchar(RNA), 3))
	peptide <- TranslateCodon(codons)
	n <- length(peptide)

	# check for stop codon
	if (identical("*", peptide[n])) {
		peptide <- peptide[1:(n - 1)]
	}

	if (!vector) {
		peptide <- paste(peptide, sep = "", collapse = "")
	}
	return(peptide)
}
TranscribeDNA("GTGAAACTTTTTCCTTGGTTTAATCAATAT")
TranscribeDNA <- function(DNA){
	RNA <- gsub("T", "U", DNA, fixed = TRUE)
	return(RNA)
}
FindPeptide <- function(DNA, peptide){
	mRNA_length <- 3 * nchar(peptide)
	encoding_sequences <- character(2 * nchar(DNA))

	for (i in 1:(nchar(DNA) - mRNA_length + 1)) {
		reading_frame <- substring(DNA, i, i + mRNA_length - 1)
		mRNA <- TranscribeDNA(reading_frame)
		pep <- TranslateRNA(mRNA)

		rev_mRNA <- TranscribeDNA(ReverseComplement(reading_frame))
		rev_pep <- TranslateRNA(rev_mRNA)

		if (identical(pep, peptide) || identical(rev_pep, peptide)) {
			encoding_sequences[i] <- reading_frame
		}
	}

	encoding_sequences <- encoding_sequences[nchar(encoding_sequences) != 0]
	return(encoding_sequences)
}
