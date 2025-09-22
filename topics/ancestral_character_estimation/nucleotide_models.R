rm(list=ls())

library(ape)
library(expm) # matrix exponential

# This code simulates output from the ace function. We
# will first try to grasp how the model behaves before we
# we try fitting it to data. 

# Ancestral character estimation methods can be based on a 
# variety of approaches, but we will consider Markov chains
# superimposed on phylogenies. In the accompanying code to 
# this we will learn how to fit these models using Maximum
# Likelihood, but for now we focus on simulating Markov chains
# on trees to understand how the models work.

# 1. Basic Markov chains

# Let's begin by considering a simple 2-state Markov chain
# defined in continuous time.

# The vector P has entries that sum to one, and the components
# are the probability the state resides in one state or the
# other at a particular time.
# Proabilities evolve over time according to the Kolmogorov equation

# dP_1/dt = - q12 P1 + q21 P2,
# dP_2/dt = q12P1 - q21 P2 

# If we treat P as a column vector,
# this can be written in vector-matrix form
# dP/dt = Q P,
# where Q is a stochastic rate matrix with columns summing to zero

# In some cases, the row vector P^T is used instead. In this case,
# the equation takes the form dP^T/dt = P^T Q

# We can solve these differential equations by hand, using the fact
# that P1 + P2 = 1

# dP1/dt = -q12(P1) + q21 * (1-P1)
# 	 = - (q12 + q21) * P1 + q21
# rearrange to:
# dP1/dt + (q12+q21) * P1 = q21

# A standard approach to solving this type of equation is to use an
# integrating factor e^{(q12+q21)*t}

# d/dt (P1 * e^{(q12+q21)*t}) = q21 * e^{(q12+q21)*t} 
# We can integrate this up once:
# P1 * e^{(q12+q21)*t} = \frac{q21}{q12+q21} * e^{(q12+q21)*t} + C

# We can now isolate P1:

# P1 =  \frac{q21}{q12+q21}  + C e^{-(q12+q21)*t} 

# We use initial conditions to determine C. Set P10 = P1(0):
# P10 = \frac{q21}{q12+q21}  + C e^{-(q12+q21)*0} 
# P10 = \frac{q21}{q12+q21}  + C 

# So now we know that
# P1 =  \frac{q21}{q12+q21}  + (P10 - \frac{q21}{q12+q21}) * e^{-(q12+q21)*t} 
# P1 = P10 * e^{-(q12+q21)*t} + \frac{q21}{q12+q21} * (1-e^{-(q12+q21)*t})
# Notice that this has a term corresponding to the initial condition that 
# is dominant for t \approx 0, and in the long run P1 \to \frac{q21}{q12+q21}

# We can also find P2, using the fact that P2 = 1-P1:
# P2 = \frac{q12+q21}{q12+q21} - \frac{q21}{q12+q21} + \frac{q21}{q12+q21}*e^{-(q12+q21)*t} - P10 *e^{-(q12+q21)*t}

# P2 = \frac{q12}{q12+q21} + (\frac{q21}{q12+q21} - P10)*e^{-(q12+q21)*t}
# P2 = \frac{q12}{q12+q21} + (\frac{q21}{q12+q21} - 1 + P20)*e^{-(q12+q21)*t}
# P2 = \frac{q12}{q12+q21} + (\frac{-q12}{q12+q21} + P20)*e^{-(q12+q21)*t}
# P2 = \frac{q12}{q12+q21}*(1-e^{-(q12+q21)*t}) + P20*e^{-(q12+q21)*t}

# So now, we have our solutions for the Markov chain probabilities over time:

P <- function(t, P0, q12,q21){
	P1 <- function(t){
		P0[1]*exp(-(q12+q21)*t) + q21/(q12+q21) * (1-exp(-(q12+q21)*t))
	}
	P2 <- function(t){
		P0[2]*exp(-(q12+q21)*t) + q12/(q12+q21) * (1-exp(-(q12+q21)*t))
	}
	P <- c(P1(t), P2(t))
	return(P)
}

# Exercises:
# 1. Write down a 4-state Markov model corresponding to A, C, T, and G 
#	in DNA. Evaluate the transition probabilities as above (hint:
#	the matrix exponential can be calculated in R using 
#	the expm function from the package of the same name

# 2. The Jukes-Cantor model is the simplest nucleotide substitution model.
#	By assumption, the equilibrium probabilities of each nucleotide
#	are equal, and the transition rates between states are likewise
#	all equal. What is the generator for this process? Write a function
#	to describe transition probabilities as a function of time. 
#	(Hint: how many parameters will your generator have
#	if all transitions between states occur at equal rates?)

# 3. The K80 (Kimura, 1980) model elaborates on the Jukes-Cantor model by
#	allowing substitutions within pyrimidines (C,T) and within purines 
#	(A, G) to occur at a different rate than "transversions" between
#	 pyrimidines and purines. In biology, substitutions between two 
#	purines or between two pyrimidines are called "transitions" (but
#	in Markov chains, mathematicians refer to all state changes as
#	transitions -- confusing). Assuming that transitions (in the 
#	biological sense of the word) occur at the same rate regardless of
#	whether the nucleotides are purines or pyrimidines, and assuming
#	that transversions in either direction occur at the same rate, 
#	write the generator for the Markov process and write a function that
#	calculates transition probabilities as a function of time.


