library(Matrix)
library(saasi)
library(diversitree)
library(tidytree)
library(ggtree)
library(ggplot2)
#library(ggimage)
library(phangorn)

set.seed(2)

# 2 trees: one ultrametric, one not

# Number of states?
nstates <- 5

pars <- data.frame(
	state=c(1:nstates),
	prior=rep(1,nstates)/nstates,
	lambda=runif(nstates,min=3,max=6),
	mu=rep(1,nstates),
	psi=c(runif(nstates, min=0.5,max=2)))

pars$lambda[1] <- 3

#q_matrix <- matrix(0.15, 2,2); diag(q_matrix) <- NA

U <- matrix(0, nstates, nstates)
U[upper.tri(U, diag = TRUE)] <- 1:nstates*(nstates+1)/2   # just an example filling

diag(U) <- 0

U[1,c(nstates-1,nstates)] <- 0
U[2,c(nstates)] <- 0

U <- U + t(U)

U <- U/50

diag(U) <- -colSums(U)

q_matrix <- U


#set.seed(1)

# create the tree object

phy <- sim_bds_tree(pars, q_matrix, x0=1, max_taxa = 300, max_t = 300,
             include_extinct = FALSE)

h <- history.from.sim.discrete(phy, 1:nstates)

true_phy_info <- as_tibble(phy)
true_phy_info$State <- c(factor(h$tip.state),factor(h$node.state))

node_depths <- node.depth.edgelength(phy)

tmrca <- max(node_depths)

# check which tips are at the present day
tips_to_drop <- phy$tip.label[abs(node_depths[1:length(phy$tip.label)] - tmrca) <= 0.01]

phypast <- drop.tip(phy, tips_to_drop)

true_phypast_info_new <- as_tibble(phypast) %>% mutate(State = c(factor(h$tip.state),factor(h$node.state))[label])
phypast$tip.state <- phypast$tip.state[setdiff(names(phypast$tip.state), tips_to_drop)]


# create the new tree by dropping the tips at the present day
phypresent <- keep.tip(phy, tips_to_drop)

true_phypresent_info_new <- as_tibble(phypresent) %>% mutate(State = c(factor(h$tip.state),factor(h$node.state))[label])
phypresent$tip.state <- phypresent$tip.state[setdiff(names(phypresent$tip.state), tips_to_drop)]


p0 <- ggtree(phy) %<+% true_phy_info + geom_point(aes(color=State),size=2) +
  ggtitle("True Phylogeny") +
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
p0


p1 <- ggtree(phypresent) %<+% true_phypresent_info_new + geom_point(aes(color=State),size=2) +
  ggtitle("True Phylogeny") +
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
p1


p2 <- ggtree(phypast) %<+% true_phypast_info_new + geom_point(aes(color=State),size=5) +
  ggtitle("True Phylogeny - without present day tips") +
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))

p2

phyaln <- simSeq(phy, l=1000, type='DNA', rate=1/100)

# uncomment to save stuff:
#save(file='phys.Rdata', phy, phypast, phypresent)
#write.phyDat(phyaln, file="phyaln.fasta", format="fasta")

# run ace using the joint and marginal reconstruction:


ace_phy <- phypast
ace_phy$node.label <- NULL

# Note: Do not have this problem if use earlier version `ape`
# Error in names(obj$ace) <- phy$node.label :
# attempt to set an attribute on NULL

#tip data:
ace_phy$tip.state <- ace_phy$tip.state[setdiff(names(ace_phy$tip.state), tips_to_drop)]

#fit ace using joint reconstruction:
asr_ace_joint<-ace(ace_phy$tip.state, ace_phy,type = "discrete", model="ER")

#extract node probabilities:
ace_node_lik <- as.data.frame(asr_ace_joint$lik.anc)
ace_node_lik$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
ace_pie <- nodepie(ace_node_lik,cols=1:5)

par(mfrow=c(1,2))
plot(phypast, show.tip.label=F)
tipstates <- phypast$tip.state
names(tipstates) <- phypast$tip.label
tiplabels(pch=21, bg=c('blue','orange','red','purple','darkgreen')[phy$tip.state[phy$tip.label]],
	cex=2)
nodelabels(pch=21, bg=c('blue','orange','red','purple','darkgreen')[phy$node.state[phy$node.label]],
	cex=2)

plot(phypast, show.tip.label=F)
tipstates <- phypast$tip.state
names(tipstates) <- phypast$tip.label
tiplabels(pch=21, 
	bg=c('blue','orange','red','purple','darkgreen')[phy$tip.state[phy$tip.label]],
	cex=2)
nodelabels(pch=21, pie=ace_node_lik, bg = c('blue','orange','red','purple','darkgreen'),
	cex=.2)


p3 <- ggtree(ace_phy) %<+% true_phypast_info_new + geom_tippoint(aes(color=State),size=2)+
  ggtitle("ace") +
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
p3 <- inset(p3, ace_pie,width = 0.07,height = 0.07,hjust=0.005)


p3

asr_ace_marginal<-ace(ace_phy$tip.state, ace_phy,type = "discrete", model="ER", marginal=TRUE)
ace_node_lik <- as.data.frame(asr_ace_marginal$lik.anc)
ace_node_lik$node <- 1:ace_phy$Nnode + Ntip(ace_phy)
ace_pie <- nodepie(ace_node_lik,cols=1:5)
p4 <- ggtree(ace_phy) %<+% true_phypast_info_new + geom_tippoint(aes(color=State),size=2)+
  ggtitle("ace") +
  theme(text = element_text(size = 15,family = "serif"),plot.title = element_text(size=15))
p4 <- inset(p4, ace_pie,width = 0.07,height = 0.07,hjust=0.005)
p4



