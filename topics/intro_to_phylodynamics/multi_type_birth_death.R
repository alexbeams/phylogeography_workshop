library(TreeSim)

#We will examine the multi-type birth death model without migration.
#Individuals in location i can give birth to individuals in any state j,
#but lineages do not move between locations otherwise.

n<-200
lambda <- rbind(c(.05,1),c(1,.2))
mu <- c(1,1)
sampprob <-c(0.1,0.1)
numbsim<-1
trees<-lapply(rep(n,numbsim),sim.bdtypes.stt.taxa,
	lambdavector=lambda,deathvector=mu,sampprobvector=sampprob)

tree <- trees[[1]]
plot(tree, show.tip.label=F)
tiplabels(pch=21,bg=c('blue','red')[tree$states])

## Exercises
## 1.
#Suppose that we want to use this as a model for phylogeography.
#The sim.bdtypes.stt.taxa function does not model migrations per se,
#but rather births from state i to state j. If migration events are
#rare relative to replication, then maybe this is a reasonable framework.
#Use ace(tree$state,tree) with any constraints on the Markov process's 
#Q matrix that you like, and plot the reconstructed ancestral states
#on the tree (unfortuantely, sim.bdtypes.stt.taxa does not give 
#information about states along edges, or at internal nodes).

## 2. 
#If the matrix $\Lambda$ has large off-diagonal terms, what do you see
#in your simulated trees? Could this be used to model anything biologically
#interesting?

## 3.
#If the matrix $\Lambda$ is diagonal-dominant but one of the $\lambda_{ii}$
#is much larger than the rest, what do you notice about the simulated 
#trees? Comment on the ability to estimate parameters of models in this case. 
