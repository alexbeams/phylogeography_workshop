rm(list=ls())

# Now that we have some familiarity with the behavior of the birth-death
#	model, we will try to use the model for parameter inference,
#	i.e. model fitting.

# Model fitting just consists of selecting parameters of the model so
#	that the simulated outputs match our data as closely as possible.
#	The dumbest way imaginable is to use the "eyeball test". The most
#	sophisticated approach is to use likelihoods. Methods like synthetic
#	likelihood and approximate Bayesian computations (ABC) fall somewhere
#	in between.

# For this set of examples and exercises, we will just adopt a frequentist 
#	approach to statistical inference.

# The diversitree package has a convenient class of functions that can be
#	used to define likelihood functions for different models. We will
#	start with the simplest case - the birth death model.

# Example 1: Benchmarking from simulations

# simulate a list of bd trees here

tree <- tree.bd(c(lambda=5,mu=0), max.taxa = 20)

loglik <- make.bd(tree)

# the make.bd function is a scalar function of two variables
#	of the form loglik = loglik(c(lambda,mu))

# Let's visualize the likelihood surface one parameter at a time 
#	(holding the other parameter at the true value):

loglik_lambda <- function(lambda) loglik(c(lambda,0))
loglik_mu <- function(mu) loglik(c(1,mu))


lambdavals <- exp( seq(-3,3,length=30))
muvals <- exp( seq(-6,2,length=30))
loglik_lambda_vals <- sapply(lambdavals, loglik_lambda)
loglik_mu_vals <- sapply(muvals, loglik_mu)

loglik_max <- max(c(loglik_lambda_vals, loglik_mu_vals))

par(mfrow=c(2,1))
plot(lambdavals, loglik_lambda_vals, ylim=c(loglik_max-10,loglik_max+2))
plot(muvals, loglik_mu_vals, ylim=c(loglik_max-10,loglik_max+2))

# The maximum of the log-likelihood function roughly coincides with the
#	parameters we used. How variable is this? Let's simulate more 
#	trees:

treelist <- trees(c(lambda=5,mu=0), 
		type=c("bd"), n=100, 
		max.taxa=20)

# make a log-likelihood function for each tree:
logliklist <- lapply(treelist, make.bd )

loglik_lambda_list <- lapply(logliklist, function(x){
			function(lambda) x(c(lambda,0))
	})

loglik_mu_list <- lapply(logliklist, function(x){
			function(mu) x(c(mu,0))
	})


lambdavals <- exp( seq(-2,2,length=30))
muvals <- exp( seq(-6,2,length=30))

loglik_lambda_list_vals <- sapply(loglik_lambda_list,
	function(x){sapply(lambdavals, x)})

loglik_mu_list_vals <- sapply(loglik_mu_list,
	function(x){sapply(muvals, x)})

loglik_max <- max(c(loglik_lambda_list_vals,loglik_mu_list_vals))
loglik_min <- min(c(loglik_lambda_list_vals,loglik_mu_list_vals))


par(mfrow=c(2,1))
plot(lambdavals, loglik_lambda_list_vals[,1]/abs(max(loglik_lambda_list_vals[,1])))
for(i in 2:100){
	lines(lambdavals, loglik_lambda_list_vals[,i]/abs(max(loglik_lambda_list_vals[,i])))
}

plot(muvals, loglik_mu_list_vals[,1]/abs(max(loglik_mu_list_vals[,1])))
for(i in 2:100){
	lines(muvals, loglik_mu_list_vals[,i]/abs(max(loglik_mu_list_vals[,i])))
}

# So, most of these have a maximum in the vicinity of the true parameter values that we
#	suppied, but not all.

# The "Maximum Likelihood Estimates" are the values of (lambda,mu) where the likelihood function is
#	maximized. What does the distribution of MLEs look like for this set of trees?

mles <- lapply(logliklist, function(x){find.mle(x, c(1,0), method='subplex')})

mles_pars <- t(sapply(mles, function(x) x$par))

par(mfrow=c(1,2))
hist(mles_pars[,1],main='',
	xlab=bquote(hat(lambda)))
hist(mles_pars[,2],main='',
	xlab=bquote(hat(mu)))

# Remember that the true value of mu we used was mu=0. A lot of the models suggest values of 
#	mu > 0, which is interesting. This tells us that, in practice, we will probably not
#	be able to reliably estimate mu when it is close to zero but positive.

# We can frame a hypothesis test. Null hypothesis: mu=0. Alternative hypothesis: mu > 0. 
#	Supposing we wish to keep false-positives at 5% or less, we could propose the statistical
#	test that we reject the null hypothesis if our estimate of mu is at the 95-th percentile or greater
#	of this "null" distribution. In our case, this critical value of mu would be approximately

mu_critical <- quantile(mles_pars[,'mu'], 0.95) 

# So in the future, if we estimate a value of mu greater than mu_critical, and we say "this appears 
#	significantly different from zero, therefore we reject the null hypothesis", we will be making a 
#	type 1 error ("false-positive") less than 5% of the time. A 1/20 chance of a mistake is reasonable
#	in some circumstances, but less so in others.

##############
# Exercises:
##############

# 1. Suppose we want a more stringent statistical test for the previous example. What is the critical
#	value of mu corresponding to a 1% chance of comitting a Type 1 error?

# 2. Modify the previous example to see how things change for a different value of lambda. Try lambda = 5.
#	Does the critical value of mu change? 

# 3. The Generalized Likelihood Ratio Test: another way to assess whether mu is significantly different from 
#	zero is to fit two versions of the log-Likelihood function: the one we already did, and a constrained
#	version of the model with the constraint mu=0 hard-coded in (so that we just estimate lambda). The 
#	model with more parameters will always fit the data better, but it may not fit the data "that much" 
#	better. The GLT looks a log(ratio of likelihoods) = difference(log-Likelihoods) to produce a test
#	statistic (also, rather annoyingly in this case, called lambda). Use constrain(loglik, mu ~ 0) to 
#	create the constrained log-likelihood functions, and fit them to obtain new estimates for lambda. 
#	Use Wikipedia, Google, or your favorite (work-safe) AI companion to help you formulate GLT test 
#	statistics and assess significant departues from the null hypothesis mu=0. Do the results of this
#	statistical test match with those from earlier? If not, why not?  
#	(Hint: the GLR also makes an asymptotic approximation. What is that approximation, and under
#	what situations would it apply here?)
 
# 4. Consider a different null hypothesis: lambda > mu. Develop a test statistic for this hypothesis, and 
#	identify critical values of your test statistic for 95% and 99% confidence regions. (Simulate large
#	numbers of trees to obtain MLEs; don't worry about GLR tests unless you want to).

# 5. If your simulation produced trees with estimated values of lambda < mu, plot these trees and comment
#	on any features they exhibit that stand out to you. If your simulation did not produce trees with
#	MLEs having lambda < mu, simulate as necessary until you find some.

# 6. The find.mle function uses a variety of numerical optimization algorithms to find MLEs. There is also a
#	function called "mcmc". Use this function with nsteps=1000 and w=c(2,1). Plot loglikelihood ("p")
#	against parameters to visualize profile log-likelihoods, and also plot lambda against mu. Comment 
#	on what you see, and discuss implications of correlations between lambda and mu. (Note: there are
#	other mcmc routines besides this one in diversitree. If you have a favorit method, feel free to 
#	use it.) 

# 7. The previous exercises and examples used the max.taxa argument to condition on trees of a particular
#	size. Flip a coin ( rbinom(1,1,0.5)). If it comes up heads (=1), change max.taxa to a different number
#	and repeat Exercises 1-6. If it comes up tails (=0), replace max.taxa with max.t to simulate trees 
#	for a fixed amount of time (with variable numbers of taxa), and then repeat Exercises 1-6. (Caution:
#	injudicious choices of lambda and max.t can make your computer explode. Ctrl + c is useful (on Mac, 
#	anyway) for cancelling calculations in the Terminal.) 


