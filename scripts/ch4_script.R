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
