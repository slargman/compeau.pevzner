# 4A
RNA <- readLines("../data/rosalind_ba4a.txt")
peptide <- TranslateRNA(RNA)
writeClipboard(peptide)

# 4B
dataset <- readLines("../data/rosalind_ba4b.txt")
DNA <- dataset[1]
peptide <- dataset[2]
encoding_seq <- FindPeptide(DNA, peptide)
writeClipboard(encoding_seq)

# 4C
peptide <- readLines("../data/rosalind_ba4c.txt")
cyclospectrum <- PeptideSpectrum(peptide)
output <- paste(cyclospectrum, sep = "", collapse = " ")
writeClipboard(output)
