library(coda) # for handling output from MCMC

# This code is complementary to the other one in this subdirectory, 
# and walks through how to fit models using a character mapping approach.

# Most of the phytools stuff is based on the Mk model, which is a particular
# type of Markov model. M stands for "Markov", and "k" indicates the number of 
# states. (Lewis has an interesting discussion relating to the use of conditional
# likelihoods as a means of overcoming acquisition bias -- the tendency to only
# record changes of state, and not maintenance of a character).


# fitMk fits the model to the data using nlminb, and optimization routine similar to
# quasi-Newton methods like BFGS. Compare this to ace.

# mcmcMk fits Mk models to trees using MCMC.


data(sunfish.tree)
data(sunfish.data)

## extract discrete character (feeding mode)
fmode<-setNames(sunfish.data$feeding.mode,
    rownames(sunfish.data))

## fit "ER" model
fit.ER<-fitMk(sunfish.tree,fmode,model="ER")
print(fit.ER)

# this will run for a while, and you should see traceplots updating
# in real time until it finishes:
fit.mcmc <- mcmcMk(sunfish.tree, fmode, model="ER", ngen=10000)

# the dimensions of fit.mcmc will be ngen * (number of parameters in the Mk model + 2),
# which for this run will be 10000 x 3

# The first column of fit.mcmc is timesteps of the MCMC routine
# The last column of fit.mcmc is the logLikelihood (logLik)
# the columns between first and last are the parameters at each timestep

# We can look at the traceplots:
plot(fit.mcmc)

# We can also plot the posterior distribution for our parameters:
# By default, this should toss out the first 20% of the MCMC updates
# as part of a burn-in (we discard transients because we only care about
# the stationary distribution of the MCMC run)
plot(density(fit.mcmc))

# Of note, we can also use base R plotting to accomplish these things as well. Here
# is a histogram of the parameter values of the entire MCMC run:
hist(fit.mcmc[,'[pisc,non]'],breaks=100) 

# and the last 80% of the run:
hist(fit.mcmc[c(2001:10000),'[pisc,non]'],breaks=100) 

## I also really like to plot profile log-likelihood (logLik~parameters):

plot(fit.mcmc[,'logLik']~fit.mcmc[,'[pisc,non]'])
## ^^ This helps us visualize the peak of the likelihood surface. Can be useful for
## diagnosing convergence problems


# We could fit a more complicated model with two rates:
fit.mcmc <- mcmcMk(sunfish.tree, fmode, model="ARD", ngen=10000)

# Now, we have profile log-likelihoods in two parameters:

par(mfrow=c(2,1))
plot(fit.mcmc[,'logLik']~fit.mcmc[,'[pisc,non]'])
plot(fit.mcmc[,'logLik']~fit.mcmc[,'[non,pisc]'])
# ^^ this one looked weird. Let's discard MCMC transients:
 
par(mfrow=c(2,1))
plot(fit.mcmc[9000:10000,'logLik']~fit.mcmc[9000:10000,'[pisc,non]'])
plot(fit.mcmc[9000:10000,'logLik']~fit.mcmc[9000:10000,'[non,pisc]'])


# We can also visualize correlations in parameters:
plot(fit.mcmc[,'[pisc,non]'], fit.mcmc[,'[non,pisc]'])

## 1. Comment on any interesting things you see in profile loglikelihoods for models with
##	many parameters. How many local optima does it look like your MCMC runs find?
## 2. 


