#read in the hMPV phylogeny and metadata tsv from Nextstrain

library(readr)
library(ape)
library(ggtree)
library(phytools)

tree <- read.tree('nextstrain_hmpv_all_genome_timetree.nwk')

dat <- read_tsv('nextstrain_hmpv_all_genome_metadata.tsv')

#tip_data <- dat$region
tip_data <- dat$clades
names(tip_data) <- dat$strain

cols <- setNames(rainbow(length(unique(tip_data))), unique(tip_data))

plotTree(tree, ftype="off")  # ftype="off" suppresses tip labels
tiplabels(pch=19, col=cols[tip_data], cex=1.1)

