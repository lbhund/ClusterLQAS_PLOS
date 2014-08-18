
######################################################
# IMPORTANT: first, run entire pezzoli.R file in R,  #
# to read the functions into R;			     #
# this must be done every time you start up R	     #
######################################################

# design survey for 75-90% couplet with risks .1, sd .1, 
# and fixing the number of sampled clusters to 6
design <- lqaspezzoli(pl=.75, pu=.9, sd = .1, alpha = 0.1, beta = 0.1, k=6, nsim=5000) 
summary(design)

# determine errors alpha and beta for 6x10 design with d=50 for 75-90% couplet 
eval <- lqasriskpezzoli(m=10, k=6, d=50, pl=.75, pu=.9, sd=.1, nsim=5000)
summary(eval)

# simulate 1000 draws from binomial-scaled distribution with mean .75, sd .1
# and plot histogram
pj <- pezzoli(p = .75, sd = .1, k = 1000)
hist(pj, main = "Histogram of cluster-level coverages")