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

# 4D
m <- as.integer(readLines("../data/rosalind_ba4d.txt"))
num <- CountPeptidesWithMass(m)
writeClipboard(as.character(num))

# 4E
dataset <- readLines("../data/rosalind_ba4e.txt")
spectrum <- as.numeric(strsplit(dataset, split = " ", fixed = TRUE)[[1]])
peptides <- CyclopeptideSequencing(spectrum)
peptides_num <- paste(peptides, sep = "", collapse = " ")
writeClipboard(peptides_num)

# 4F
dataset <- readLines("../data/rosalind_ba4f.txt")
peptide <- dataset[1]
spectrum <- as.numeric(strsplit(dataset[2], split = " ", fixed = TRUE)[[1]])
score <- as.character(PeptideScore(peptide, spectrum))
writeClipboard(score)

# 4G
dataset <- readLines("../data/rosalind_ba4g.txt")
N <- as.integer(dataset[1])
spectrum <- as.numeric(strsplit(dataset[2], split = " ", fixed = TRUE)[[1]])
peptide <- LeaderboardCyclopeptideSequencing(spectrum, N)
writeClipboard(peptide[1])

# 4H
dataset <- readLines("../data/rosalind_ba4h.txt")
spectrum <- as.numeric(strsplit(dataset[1], split = " ", fixed = TRUE)[[1]])
convolution <- SpectrumConvolution(spectrum)
convolution <- paste(convolution, sep = "", collapse = " ")
writeClipboard(convolution)
