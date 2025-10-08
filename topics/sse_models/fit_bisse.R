rm(list=ls())

library(diversitree)
library(mcmcensemble)


# Now that we have some familiarity with the behavior of the BISSE
#	model, we will try to use the model for parameter estimation,
#	i.e. inference via model fitting.

# Goals:
#	1. Learn how to specify a BiSSE model log-likelihood
#	2. Learn how to fit a BiSSE model
#	3. Learn how to constrain models to estimate fewer parameters
#	4. Develop an intuition for how well we can estimate parameters



# Example: Benchmarking from simulations

# simulate a BiSSE tree here

pars <- c(
        lambda0=1,
        lambda1=5,
        mu0=0.8,
        mu1=0.8,
        q01=0.1,
        q10=0.12)

tree <- tree.bisse(pars, max.taxa=50, x0=0)

# the make.bisse function is a scalar function of two variables
#	of the form loglik = loglik(c(lambda,mu))

loglik <- make.bisse(tree, tree$tip.state)

# What are the MLEs?
mles <- find.mle(loglik, pars)

# these are the parmeter values that this infers:
mlepars <- mles$par

# Let's visualize the likelihood surface one parameter at a time
#	(holding the other parameter at the MLE; not exactly profile
#	likelihoods, but pretty close):

loglik_lambda0 <- function(lambda0) loglik(c(lambda0,mlepars[-1]))
loglik_lambda1 <- function(lambda1) loglik(c(mlepars[1],lambda1,mlepars[3:6]))

loglik_mu0 <- function(mu0) loglik(c(mlepars[1:2],mu0,mlepars[4:6]))
loglik_mu1 <- function(mu1) loglik(c(mlepars[1:3],mu1,mlepars[5:6]))

loglik_q01 <- function(q01) loglik(c(mlepars[1:4],q01,mlepars[6]))
loglik_q10 <- function(q10) loglik(c(mlepars[1:5],q10))

lambdavals <- exp( seq(-5,2,length=50))
muvals <- exp( seq(-9,2,length=50))
qvals <- exp( seq(-9,2, length=50) )

loglik_lambda0_vals <- sapply(lambdavals, loglik_lambda0)
loglik_lambda1_vals <- sapply(lambdavals, loglik_lambda1)

loglik_mu0_vals <- sapply(muvals, loglik_mu0)
loglik_mu1_vals <- sapply(muvals, loglik_mu1)

loglik_q01_vals <- sapply(qvals, loglik_q01)
loglik_q10_vals <- sapply(qvals, loglik_q10)

par(mfrow=c(3,2))

plot(lambdavals,loglik_lambda0_vals)
abline(v=pars['lambda0'],col='red')

plot(lambdavals,loglik_lambda1_vals)
abline(v=pars['lambda1'],col='red')

plot(muvals,loglik_mu0_vals)
abline(v=pars['mu0'],col='red')

plot(muvals,loglik_mu1_vals)
abline(v=pars['mu1'],col='red')

plot(qvals,loglik_q01_vals)
abline(v=pars['q01'],col='red')
plot(qvals,loglik_q10_vals)
abline(v=pars['q10'],col='red')

# So, varying each parameter while holding the others at the values used
#	to produce the simulated BiSSE tree indicate MLEs may be a good
#	approach to estimating parameters.

# The MLE differs from the "true" values used to produce the simulation, but
#	that is expected.

# Of course, the profile likelihoods don't show how correlated some of these
#	parameters are with each other. We might as well use MCMC to do this:

# This will take quite a lot longer than our birth-death model with
#	just 2 parameters:
#mcmc_fit <- mcmc(loglik, pars, nsteps = 1000, w=1)


# let's use an ensemble MCMC sampler:
# make.bisse will throw errors if we try to plug in negative parameters:
#loglik(-pars)

# For an MCMC routine that isn't included in diversitree, we need to ensure that
#	we can do an unrestricted search in parameter space.

loglik_transformed <- function(pars){
	pars <- exp(pars)
	loglik <- loglik(pars)
	return(loglik)
}

loglik_transformed(log(pars))
# This looks ok.

loglik_transformed(-log(pars))
# Calculates a value without error. Good to go.

nwalkers <- 100 #ensemble size
nsteps <- 50 #number of times to update entire nsemble

# initialize the ensemble:
mcmc_inits <- matrix(nrow=nwalkers, ncol=6)
for(i in 1:nwalkers) mcmc_inits[i,] <- jitter(log(mlepars))

# run the ensemble MCMC:
mcmc_fit <- MCMCEnsemble(loglik_transformed,
	inits=mcmc_inits,
	max.iter=nwalkers*nsteps,
	n.walkers = nwalkers)

# Now, let's try to visualize the profile log-likelihoods (including
#	the curves we generated earlier):

par(mfrow=c(3,2))

plot(mcmc_fit$log.p ~ mcmc_fit$samples[,,1],xlab=bquote(lambda[0]))
#lines(log(lambdavals), loglik_lambda0_vals, col='red')

plot(mcmc_fit$log.p ~ mcmc_fit$samples[,,2],xlab=bquote(lambda[1]))
#lines(log(lambdavals), loglik_lambda1_vals, col='red')

plot(mcmc_fit$log.p ~ mcmc_fit$samples[,,3],xlab=bquote(mu[0]))
#lines(log(muvals), loglik_mu0_vals, col='red')

plot(mcmc_fit$log.p ~ mcmc_fit$samples[,,4],xlab=bquote(mu[1]))
#lines(log(muvals), loglik_mu1_vals, col='red')

plot(mcmc_fit$log.p ~ mcmc_fit$samples[,,5],xlab=bquote(q[0~1]))
#lines(log(qvals), loglik_q01_vals, col='red')

plot(mcmc_fit$log.p ~ mcmc_fit$samples[,,6],xlab=bquote(q[1~0]))
#lines(log(qvals), loglik_q10_vals, col='red')

# In the vicinity of the MLE, the MCMC results should closely match
#	the profile likelihood curves (but far away from the MLE there
#	is no reason why they should be similar). (Why?)

# Looking at the different plots, are you able to say which parameters
#	that are more difficult to estimate than others?


# What if we look at parameter correlations?
#plot(mcmc_fit$samples[,,1] ~ mcmc_fit$samples[,,2],
#	xlab=bquote(lambda[1]),
#	ylab=bquote(lambda[0]))

# Exercise 1. What parameters seem to exhibit the highest correlations? The lowest?
#	(You will need to plot all the pairwise combinations.) How do parameter
#	estimates compare to the values used in the simulation? Is it easier to
#	estimate one diversification rate or the other? Ditto for the mu's and q's.
#	Why do you think that is?

# Exercise 2. Increase the number of taxa in your tree. How does your ability
#	to estimate parameters change? Does it become easier to estimate migration
#	rates, or death rates? Using 200 taxa should be doable, but see how high
#	you can go before computations get bogged down. Be sure to describe
#	bias as well as precision of estimates.

# Exercise 3. You may notice that the "data" (putting it in quotes b/c we simulated it)
#	contain more information about one of the migration rates than the other. Which
#	migration rate is it, and why do you think it might be easier to estimate
#	it with greater precision than the other one? It might be helpful to plot
#	the tree.

#plot(history.from.sim.discrete(tree, states=c(0,1)),
#                tree, col=c('0'='black','1'='red'))

# Exercise 4. Set migrations to be the same order of magnitude as the
#	the diversification rates and repat the analysis. Does it always
#	become easier to estimate q_ij?

# Exercise 5. Set migration to be asymmetric. How easy is this to detect?
#	You can start by assuming diversification and extinction rates
#	are the same in each location.

# Exercise 6. Under what circumstances do you think it is possible to reliably
#	estimate extinction rates? Carry out an analysis to confirm or refute your
#	hypothesis.

# Exercise 7. Are there are any situations where it is difficult to estimate
#	diversification rates? Comment on bias and variance.

# Exercise 8. Naively, one might think that having larger differences between
#	parameters (lambda0 vs lambda1, for example) would make estimation easier.
#	Why is that not necessarily the case with BiSSE models?
#	Think back to aeons ago when we simulated these models under different
#	parameter combinations...

# Exercise 9. Try setting a constraint in the log-likelihood function that mu0=mu1.
#	Does this help estimate the other parameters? Does it make it easier to
#	estimate the overall exctinction rate (mu)?

# Hint: it works like this: loglik_constrained <- constrain(loglik, mu0 ~ mu1)

# Exercise 10. Load the phangorn package, and use the nni function to get all of
#	the trees that are one nni move away from your simulated tree. Evaluate
#	the likelihood of all of the trees (but just use the MLE for the other
#	parameters). Does the true tree carry the highest likelihood? Do you
#	think the tree itself is identifiable? (Hint: use the chronos function
#	to produce ultrametric trees from the nni output).
