rm(list=ls())

# By the end of this set of exercises and examples, we will code up a 
#	Metropolis-Hastings algorithm to carry out our own MCMC.

# We will start by introducing some basic Monte Carlo techniques.

# 1. Monte Carlo sampling to approximate Pi
# 2. Monte Carlo Integration 
# 3. Importance sampling
# 4. Metropolis-Hastings and MCMC

###############
# Example 1: Monte Carlo sampling to approximate Pi 
###############

# From calculus, we know how to evaluate integrals using 
# antiderivatives. You might also be familiar with numerical
# quadrature techniques like the trapezoid rule and Simpson's
# rule that approximate areas under curves with shapes that have
# nice are formulas.

# It looks like this:
# code up example trapezoid/simpsons' rule w/ picture

# Monte Carlo integration works a little bit differently. Let's start
# by approximating the area of a circle to approximate pi.

# Suppose there is a square dart board with a circle inscribed. My dart
# throws are random, so when I launch a dart at the board it lands randomly
# someplace in the square. Sometimes darts will land inside the circle, 
# sometimes they will land outside. 

# If I've chosen my square so that it is 1 meter x 1 meter, then the 
# inscribed circle will have a diamter of 1 meter. Like this:

drawDartboard <- function(){
	# Square with Inscribed Circle in R
	# Set up the plot
	plot(c(-1.2, 1.2), c(-1.2, 1.2), 
	     type = "n", 
	     asp = 1, 
	     xlab = "", 
	     ylab = "", 
	     main = "Square with Inscribed Circle")

	# Define square side length
	side_length <- 1
	half_side <- side_length / 2

	# Draw the square
	# Square vertices: bottom-left, bottom-right, top-right, top-left, back to bottom-left
	square_x <- c(-half_side, half_side, half_side, -half_side, -half_side)
	square_y <- c(-half_side, -half_side, half_side, half_side, -half_side)
	lines(square_x, square_y, col = "blue", lwd = 2)

	# Draw the inscribed circle
	# For a square with side length s, the inscribed circle has radius r = s/2
	radius <- side_length / 2
	theta <- seq(0, 2*pi, length.out = 100)
	circle_x <- radius * cos(theta)
	circle_y <- radius * sin(theta)
	lines(circle_x, circle_y, col = "red", lwd = 2)

	# Add a grid for reference
	grid(col = "lightgray", lty = "dotted")

	# Add legend
	legend("topright", 
	       legend = c("Square", "Inscribed Circle"), 
	       col = c("blue", "red"), 
	       lwd = 2)

	# Optional: Add center point
	points(0, 0, pch = 19, col = "black", cex = 0.8)
} 


# Now, I'm going to throw darts at this target. Because they land randomly,
# the (x,y) coordinates of the darts are each uniformly distributed, like this:

throwDarts <- function(ndarts){

	# randomly choose (x,y) coordinates independently:
	x <- runif(ndarts,min=-0.5, max=0.5)
	y <- runif(ndarts,min=-0.5, max=0.5)
	
	# re-draw the dartboard:	
	drawDartboard()
	
	# add the hit locations to the dartboard:
	points(x,y, pch=4)

	return(list(x=x,y=y))
} 

# Play around with the throwDarts function for different values of ndarts.
# How does this help us evaluate pi? Remember that Area(circle) = pi * radius^2,
# and Area(square) = radius^2 = 1, because we deliberately chose radius=1.

# So, the area of the circle is just pi, and the area of the rectangle is 1.

# Another way of saying this is that pi = Area(circle)/Area(rectangle).

# So what? Well, the darts give us a way to calculate the ratio of areas. 
# Convince yourself that (number of darts inside circle)/(total darts thrown)
# is equal to this ratio of areas.

# Another way of saying this is: the probability a dart lands inside the circle
# is just the ratio of Area(circle)/Area(square).

# Let's modify throwDarts to calculate the fraction of darts falling inside
# the circle:

countDartsInside <- function(ndarts){
	darts <- throwDarts(ndarts)

	x = darts$x
	y = darts$y
	r = sqrt(x^2 + y^2)
	
	numinside = sum(r < 0.5)
	return(numinside/ndarts)
} 

# The area of the circle is pi/4, the area of the square is 1

# So, 4 x (fraction of darts landing inside) = pi (approximately)


calculatePi <- function(ndarts){
	
	ratio = countDartsInside(ndarts)
	pi_approx = 4 * ratio
	
	return(pi_approx)
}


# Optional coding exercise: can you modify the code above
# to color in the darts according to whether they fall inside
# or outside the circle?

# Question: how many darts do you need to throw to approximate pi
# to three significant figures? 

# It seems like this approximation is doing.. ok. Against quadrature
# methods, this will not be nearly as computationally efficient. However,
# it is easy and straightforward to use, and generalizes well to more
# complicated tasks when quadrature methods are not tractable. 

##################
# Exercise 2: Using Monte Carlo to calculate definite integrals
##################

# The same exact trick works for calculating definite integrals, i.e.
# areas under curves. Let's calculate integrals of the following functions
# with limits of integration of a = -5 and b = 5.

#limits of integration:
a <- -5
b <- 5

# functions to integrate:
f <- function(x)sin(x)
g <- function(x) x^2
h <- function(x) 1/sqrt(2*pi) * exp(-x^2/2)

# antiderivatives of these functions (should have +C):
bigf <-function(x) -cos(x) 
G <- function(x) x^3/3
H <- function(x) pnorm(x,mean=0,sd=1)

# So, using the Fundamental Theorem of calculus, part 2, we know our
# definite integrals should evaluate to:

# bigf(5) - bigf)-5) = 0 (also b/c sin(x) is odd)
# G(5) - G(-5) = 250/3 = 83.333..
# H(5) - H(-5) = 0.9999994

# Let's plot these:
xvals <- seq(-10,10,length=1000)
par(mfrow=c(3,1))
plot(xvals, sapply(xvals, f), type='l',xlab='x',ylab='f(x) = sin(x)')
plot(xvals, sapply(xvals, g), type='l',xlab='x',ylab='g(x) = x^2')
plot(xvals, sapply(xvals, h), type='l',xlab='x',ylab='h(x) = gaussian')

# We could inscribe these shapes within rectangles (like with the dartboard 
# example), but the annoying difficulty with that is that our functions 
# have different heights. So, we would need to customize our dartboard 
# in each case. We don't want to do that. Instead, we want to modify our method
# to accomodate any function whatsoever, particularly as we might not always
# be able to "draw" the function conveniently (if it is a function of 1000's of
# variables, for example).

# Instead, we can go back to calculus. At some point in calculus, you were 
# taught "The Mean Value Theorem for Integrals". You may not remember it,
# but I 1000% guarantee that it was taught to you if you took calculus.

# It says this:

# \bar{f} = \frac{1}{b-a} \int_a^b f(x) dx.

# What this means is: each of these curves might go up or down along the
# interval from a to b, but it has an "average" height \bar{f}.

# This average is calculating by evaluating the area under the curve, 
# and then dividing that by the width of your integration interval.

# If f(x) is a constant function, then \int_a^b f dx = (b-a) * f, 
# and \bar{f} = f. Or: if you have a rectangle (the area) and divide
# that area by the width, you get the height. 

# What does this have to do with Monte Carlo integration? Here's how: we 
# are going to re-write the Mean Value Theorem for Integrals:

# \bar{f} = \frac{1}{b-a} \int_a^b f(x) dx
#		 <->
# \int_a^b f(x) dx = (b-a) \bar{f}

# If we could evaluate \bar{f} somehow, then multiplying by (b-a) would
# give us our integral. We can't evaluate \bar{f} exactly, but we can
# certainly approximate it. 

# All we have to do is this: randomly choose x-coordinates between
# a and b, evaluate f(x) at each of those points, and then average those
# values. This will give us a Monte Carlo approximation to the average 
# value of f(x):

ndarts <- 50000
xrands <- runif(ndarts, min=a, max=b)

fbar_approx <- mean( f(xrands) )
gbar_approx <- mean( g(xrands) )
hbar_approx <- mean( h(xrands) )

# Once we have our approximations to the average values of the functions
# (i.e. approximations to the heights), we just multiply those by the 
# widths on the x-axis (b-a) to get our approximation to the integral (i.e.
# the area under the curve):

integral_f <- fbar_approx * (b-a)
integral_g <- gbar_approx * (b-a)
integral_h <- hbar_approx * (b-a)

# Experiment with different values of ndarts to see how the approximations
# change. 

# Exercise:
# Modify the code above to evaluate integrals of the multivariate function
# f(x,y) = sin(x) * cos(y) 


#######################
## Importance sampling
#######################

# The mathematician in you might be annoyed at the crudeness of these 
# approximations. After all, standard numerical quadrature approaches 
# that we learn about in calculus are far more accurate. However, I hope
# you can appreciate the simplicity of this approach. Quadrature methods
# often need to be tailored for particular situations to work optimally.
# The code above can be used with very minor modification to evaluate any
# integral whatsoever.

# This does have one rather extreme drawback, however. Let's return to the 
# Gaussian example in the previous section. Suppose instead I had asked
# you to calculate the definite integral with a=5, b=infinity. 

# Let's plot this just to see what the trouble might be:

plot(xvals, sapply(xvals, h), xlab='x', ylab='Gaussian',
	type='l')
points(x=5, y=0,col='red')

# There are some challenges here. I probably can't use runif anymore 
# because it requires a maximum value that is finite. 

# Now, you might be saying to yourself, "Well, R has the ability to simulate
# variables from a standard normal distribuion, and our h(x) is just the PDF
# of a standard normal distribution. So, I could just simulate a large
# number of standard normal variates, and evaluate the fraction that
# are larger than 5. 

# "Okay," I reply. "Let's do that."

tailsample <- rnorm(1e4)

# While we're at it, let' just go ahead and use pnorm to figure out what the
# true answer really is:

tailprob = 1-pnorm(5)

print(tailprob)

# What proportion of tailsample was greater than 5?
sum(tailsample > 5)/length(tailsample)

# I got zero. Let's just use a larger sample:

tailsample <- rnorm(1e6)
sum(tailsample > 5)/length(tailsample)

# I'm seeing that I usually get 0 or 1 value larger than 5.



