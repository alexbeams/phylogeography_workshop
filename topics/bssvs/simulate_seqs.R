library(phangorn)
library(TreeSim)

seqs <- simSeq(tree, l=1600, type='DNA', rate=1/1000)

aln <- as.DNAbin(seqs)

write.phyDat(aln, file="simulated_alignment.fasta", format="fasta")


