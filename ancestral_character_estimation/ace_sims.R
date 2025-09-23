rm(list=ls())
library(ape)
library(phangorn)

# first let's assume we have a tree (a coalescent tree, in this case):
tree <- rcoal(40)

# Transition rate matrix (two states A,B)
Q <- matrix(c(-1,1,
              1,-1),2,2)
rownames(Q) <- colnames(Q) <- c("A","B")

sim = simSeq(tree, l = 1, Q=Q, bf = NULL, rootseq = NULL, 
	type = "USER", levels= c("A", "B"), ancestral=TRUE)


