# TiPS tutorial
library(TiPS)

time <- c(0, 20)

dT <- 0.001


safe_run <- function(f, ...) {
  out <- list()
  while(! length(out)) {out <- f(...)}
  out
}


# multiple patch model:

reactions <- c(
	"0 [beta1 * I1] -> I1",
	"0 [beta2 * I2] -> I2",
	"0 [beta3 * I3] -> I3",
	"0 [beta4 * I4] -> I4",
	"0 [beta5 * I5] -> I5",
	"0 [beta6 * I6] -> I6",
	"0 [beta7 * I7] -> I7",
	"0 [beta8 * I8] -> I8",
	"0 [beta9 * I9] -> I9",
	"0 [beta10 * I10] -> I10",

	"I1 [gamma1 * I1] -> 0",
	"I2 [gamma2 * I2] -> 0",
	"I3 [gamma3 * I3] -> 0",
	"I4 [gamma4 * I4] -> 0",
	"I5 [gamma5 * I5] -> 0",
	"I6 [gamma6 * I6] -> 0",
	"I7 [gamma7 * I7] -> 0",
	"I8 [gamma8 * I8] -> 0",
	"I9 [gamma9 * I9] -> 0",
	"I10 [gamma10 * I10] -> 0",

	"I1 [q21 * I1] -> I2",
	"I2 [q32 * I2] -> I3",
	"I3 [q43 * I3] -> I4",
	"I4 [q54 * I4] -> I5",
	"I5 [q65 * I5] -> I6",
	"I6 [q76 * I6] -> I7",
	"I7 [q87 * I7] -> I8",
	"I8 [q98 * I8] -> I9",
	"I9 [q109 * I9] -> I10",

	"I2 [q12 * I2] -> I1",
	"I3 [q23 * I3] -> I2",
	"I4 [q34 * I4] -> I3",
	"I5 [q45 * I5] -> I4",
	"I6 [q56 * I6] -> I5",
	"I7 [q67 * I7] -> I6",
	"I8 [q78 * I8] -> I7",
	"I9 [q89 * I9] -> I8",
	"I10 [q910 * I10] -> I9"
)

bd_simu <- build_simulator(reactions)


initialStates <- c(I1 = 1, I2 = 1, I3 = 0, I4 = 0,
	I5 = 0, I6 = 0, I7 = 0, I8 = 0, I9 = 0, I10 = 0)


time <- c(1975, 1998, 2018)


theta <- list(
	beta1 = 1,
	beta2 = 2,
	beta3 = 3,
	beta4 = 4,
	beta5 = 5,
	beta6 = 6,
	beta7 = 7,
	beta8 = 8,
	beta9 = 9,
	beta10 = 10,
	gamma1 = 1,
	gamma2 = 2,
	gamma3 = 3,
	gamma4 = 4,
	gamma5 = 5,
	gamma6 = 6,
	gamma7 = 7,
	gamma8 = 8,
	gamma9 = 9,
	gamma10 = 10,
	q12 = 1,
	q23 = 1,
	q34 = 1,
	q45 = 1,
	q56 = 1,
	q67 = 1,
	q78 = 1,
	q89 = 1,
	q910 = 1,
	q21 = 1,
	q32 = 1,
	q43 = 1,
	q54 = 1,
	q65 = 1,
	q76 = 1,
	q87 = 1,
	q98 = 1,
	q109 = 1
)

dT <- 0.01


safe_bd_simu <- function(...) safe_run(bd_simu, ...)


trajbd_tl <- safe_bd_simu(
  paramValues = theta,
  initialStates = initialStates,
  times = time,
  method = "approximate",
  tau = 0.001)

dates_bd <- seq(from=2015, to=2018, length.out=100)


bd_tree <- simulate_tree(
  simuResults = trajbd_tl,
  dates = dates_bd,
  deme = paste0('I',1:10) ,
 # sampled = c(I1 = 0.2, I2 = 0.8), # the type of individuals that are sampled and their proportion of sampling
  root = "I1", # type of individual at the root of the tree
  isFullTrajectory = FALSE, # deads do not generate leaves
  nTrials = 3,
  addInfos = TRUE) # additional info for each node

bd_tree


inode_cols <- ifelse(grepl(x=bd_tree$node.label,pattern="I2"),"blue","red")
ape::plot.phylo(bd_tree, root.edge = T, no.margin = F, align.tip.label = T)
nodelabels(pch=20,col=inode_cols)

