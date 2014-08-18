

########################################################################################################
# PARAMETERS DEFINITIONS FOR FUNCTIONS BELOW ###########################################################
# pl - lower threshold														 #
# pl - upper threshold														 #
# sd - standard deviation of cluster level coverage									 #
# alpha - risk defined in equation 1											 #
# beta - risk defined in equation 1												 #
# m - number of individuals to sample per cluster									 #
# k - number of clusters to sample (only specify one of m or k)							 #
# nsim - the code is based on simulation, and this specifies the number of sims, increase for accuracy #
########################################################################################################

#########################################
# Design survey with parameters pl, pu, #
# alpha, beta; standard deviation sd;   #
# and fixing either k clusters or m per #
# cluster sampled                       #
#########################################
lqaspezzoli <- function (pl, pu, sd, alpha = 0.1, beta = 0.1, m = NULL, 
	k = NULL, nsim = 3000, add = 1) {
    stop <- F
    fixk <- TRUE
    if (is.null(k) == T) 
        fixk <- FALSE
    if (is.null(m) == T & fixk == FALSE) 
        stop("Must specify m or k")
    if (fixk == TRUE) 
        m <- 0
    if (fixk == FALSE) 
        k <- 0
        while (stop == F) {
            if (fixk == TRUE) 
                m <- m + add
            if (fixk == FALSE) 
                k <- k + add
            n <- m * k
            d <- (round(0.5 * pl * n)):(round(2 * pu * n))

		lowervec <- sapply(1:nsim, function(a) sum(rbinom(k, m, p=pezzoli(pl, sd, k))))
		uppervec <- sapply(1:nsim, function(a) sum(rbinom(k, m, p=pezzoli(pu, sd, k))))

            alpha_vec <- sapply(d, function(dd) {
			return(length(which(uppervec <= dd))/nsim)
            })
		beta_vec <- sapply(d, function(dd) {
			return(length(which(lowervec > dd))/nsim)
            })
            ind <- rep(0, length(d))
            ind[alpha_vec <= alpha & beta_vec <= beta] <- 1
            if (sum(ind) > 0) 
                stop <- T
            rule <- d[which(ind == 1)]
	}
      out 		<- NULL
      out$rule 	<- rule
      out$n 	<- n
      out$alpha 	<- alpha_vec[which(ind == 1)]
      out$beta 	<- beta_vec[which(ind == 1)]
      out$alpha_max 	<- alpha
      out$beta_max 	<- beta
      out$pl 	<- pl
      out$pu 	<- pu
      out$m 	<- m
      out$k 	<- k
      out$sd 	<- sd
      class(out) 	<- "pezzolilqas"
      return(out)
}

summary.pezzolilqas <-function (object) {

	if(length(object$rule) > 1) 
		cat("NOTE: More than one decision rule found. \n")

    	cat("The LQAS design parameters are: \n\n")
    	cat("Sample Size:", object$n, "\n")
    	if(length(object$rule) ==1)
  		cat("Decision Rule:", object$rule, "\n\n")
    	if(length(object$rule) > 1)
		cat("Decision Rules:", object$rule, "\n\n")

    	cat("Classify as high if X > d.", "\n\n", sep = "")
    	cat("The true error levels for this design are: \n alpha =", 
        	round(object$alpha, digits = 4), "\n  beta =", round(object$beta, digits = 4), "\n\n")

    	cat("Design Parameters:", "\n")
    	cat("pl = ", object$pl, "and pu = ", object$pu, "\n")
    	cat("alpha = ", object$alpha_max, "and beta = ", object$beta_max, "\n")


	cat("Probabilities were calculated assuming cluster sampling was used.  \n")
	cat("Sample ", object$k, " clusters and ", object$m, " individuals per cluster, assuming an SD of ", object$sd, ".\n", sep="")

}


#########################################
# Calculate alpha and beta risks for a  #
# Pezzoli design with m per cluster,    #
# k clusters, rule d, thresholds pl and #
# pu, and standard deviation sd.        #
#########################################
lqasriskpezzoli <- function (m, k, d, pl, pu, sd, nsim=1000) {
    	if (pl >= pu) 
        	stop("pu must be greater than pl")

 	alpha_vec 	<- mean(sapply(1:nsim, function(ii) ifelse(sum(rbinom(k, m, p=pezzoli(pu, sd, k))) <= d, 1, 0)))
    	beta_vec 	<- mean(sapply(1:nsim, function(ii) ifelse(sum(rbinom(k, m, p=pezzoli(pl, sd, k))) >  d, 1, 0)))

  	out 		<- NULL
    	out$rule 	<- d
    	out$n 	<- m*k
	out$m 	<- m
	out$k 	<- k
    	out$pl 	<- pl
    	out$pu 	<- pu
	out$sd 	<- sd
    	out$alpha 	<- alpha_vec
    	out$beta 	<- beta_vec
    	class(out) 	<- "pezzolilqas"
    	return(out)
}


#########################################
# Simulate k cluster-level coverages    #
# from the binomial-scaled distribution #
# with mean p and standard deviation sd.#
#########################################
pezzoli 	<- function(p=.5, sd=.1, k=10){
	#p is prob threshold, sd is sd of p, k is number of clusters
	if(p*(1-p) < sd^2){
		warning("sd too large, eta set to 1.") 
	}
	n 	<- max(1, round(p*(1-p)/sd^2))
	out <- rbinom(k, n, p)/n
	return(out)
}

