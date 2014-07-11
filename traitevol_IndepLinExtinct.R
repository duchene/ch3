
require(phangorn)
require(MASS)
require(Matrix)

# The following function simulates phylogenetic structures (not yet of class phylo) according to the relationsips between a trait and diversification rate and the same trait and the probability of a substitution occurring.

# Arguments:	stepsize is the size of each time step in time (it influences the probability of an event of speciation, exinction, or substitution occurring).	branchstop is number of branches desired and stops the simulation. 	seqlen is the length of the genetic sequence to be generated.	traitstart is the initial value of the trait.			trait.r is the rate of change of the trait, which is adjusted when a trend is desired. 		FUNspr and FUNmu are the functions that define the relationship between the trait and the probability of bifurcation and substitution respectively.	Pext is the constant background probability of extinction.	D is the variance of trait evolution, which is taken as being constant.	molerror and sprerror are the error to be introduced to each of speciation and mutation probability at each step.


tr.mu.sp <- function(stepsize = 0.1, branchstop = 200, seqlen = 2000, traitstart = 50, trait.r = 0, regcoefmu = 0.01, regcoefspr = 0.02, Pext = 0.05, Dsd = 0.001, molerror = 0.001, sprerror = 0.01, direct = F, covariance = 0.0000099, meanmu = 0.009, meanspr = 0.1, q = matrix(rep(0.1, 16), 4, 4), age = 50){

# The following is a matrix where the columns are: the parent node, the daughter node, the branch length in time, and etiher the trait value or the means for mu and spr.

	if(direct == F){

		edgetable <- matrix(data = c(0, 1, 0, traitstart), nrow = 1, ncol = 4)

	} else {

		edgetable <- matrix(data = c(0, 1, 0, meanmu, meanspr), nrow = 1, ncol = 5)

	}

    substitutions <- list(matrix(c(0,0,0), 1, 3))
    
    subsPmat <- as.matrix(expm(q * stepsize)); colnames(subsPmat) <- c("a", "c", "g", "t"); rownames(subsPmat) <- c("a", "c", "g", "t");
    
    sequence <- sample(c("a", "c", "g", "t"), seqlen, replace = T)

	extinct <- 1

	time <- 1
	
	tips <- 1

	b1 <- regcoefmu

	b2 <- regcoefspr

	muvar <- molerror^2

	spvar <- sprerror^2

	 while(((tips - length(extinct)) * 2) < branchstop){

	 	    # Here we are setting constant time step size. Otherwise, use rexp().

		    dt <- stepsize
        	rt <- trait.r * dt

        	# The following restarts the simulation if the age of the tree exceeds 500.
        	time <- time + dt
        	if(time > 501){ time <- 1; extinct <- 1; substitutions <- list(matrix(c(0,0,0), 1, 3)); tips <- 1;

        		if(direct == F){

        			edgetable <- matrix(data = c(0, 1, 0, traitstart), nrow = 1, ncol = 4); print("A phylogeny grew too old!")

        		} else {

        			edgetable <- matrix(data = c(0, 1, 0, meanmu, meanspr), nrow = 1, ncol = 5); print("A phylogeny grew too old")

        		}

        	}
        	
        	tips <- length(edgetable[, 2][which(!edgetable[, 2] %in% edgetable[, 1])])

		    # If interested on the rate of evolution of the trait, we should set D differently, perhaps with dexp().

		    for(i in (1:nrow(edgetable))[if(length(extinct) > 1){ -extinct } else { 1:nrow(edgetable) }]){

		    	D <- rexp(1, Dsd)

		    	if(!(edgetable[i, 2] %in% edgetable[, 1])){ # Grab a tip lineage

					e1 <- rnorm(1, 0, molerror)

					e2 <- rnorm(1, 0, sprerror)

					FUNspr <- function(x) abs(0.1 * (1 - exp(-b2 * x)) + 0.01 + e2)

					FUNmu <- function(x) abs(0.02 * (1 - exp(-b1 * x)) + 0.001 + e1)

					if(b2 == 0){
						FUNspr <- function(x) abs(0.1 * (1 - exp(-b2 * x)) + meanspr + e2)
					}
					if(b1 == 0){
						FUNmu <- function(x) abs(0.02 * (1 - exp(-b1 * x)) + meanmu + e1)
					}

					if(direct == F){

						# The following makes sure that trait values do not become lower than roughly 2.

						if(edgetable[i, 4] < 2) edgetable[i, 4] <- 2 + rnorm(1, 0, sqrt(2 * D * dt))

		    			spr.new <- FUNspr(edgetable[i, 4])

		    			mu.new <- FUNmu(edgetable[i, 4])

		    		} else {

		    			if(edgetable[i, 4] < 0.001) edgetable[i, 4] <- 0.001 + rnorm(1, 0, 0.0001)

		    			if(edgetable[i, 5] < 0.01) edgetable[i, 5] <- 0.01 + rnorm(1, 0, 0.001)

		    			rand <- mvrnorm(1, c(edgetable[i, 4], edgetable[i, 5]), matrix(c(muvar, covariance, covariance, spvar), 2, 2))

		    			while(any(rand < 0)) rand <- mvrnorm(1, c(edgetable[i, 4], edgetable[i, 5]), matrix(c(muvar, covariance, covariance, spvar), 2, 2))

		    			spr.new <- rand[2]

		    			mu.new <- rand[1]

		    		}

		    		# The following is the total probability of an event occurring, which includes twice the equation for the rate of an event occurring.

		    		totalP <- ((spr.new + (mu.new * seqlen) + Pext) * dt) * exp(-((spr.new + (mu.new * seqlen) + Pext) * dt))

						# The following condition determines whether and event occurs.

		    		if(runif(1) <= totalP){
		    			
		    			### Define probability of events and execute them.

		    			extevent <- rbinom(1, 1, Pext)

						if(extevent == 1){

							# An extinction event occurs and a time step is added:

							extinct <- append(extinct, i)

							if(length(extinct) > length(!(edgetable[, 2] %in% edgetable[, 1]))) {

								if(direct == F){

									edgetable <- matrix(data = c(0, 1, 0, traitstart), nrow = 1, ncol = 4); extinct <- 1 ; time <- 1; substitutions <- list(matrix(c(0,0,0), 1, 3)); tips <- 1; print("A total extinction happened!")

								} else {

									edgetable <- matrix(data = c(0, 1, 0, meanmu, meanspr), nrow = 1, ncol = 5); extinct <- 1 ; time <- 1; substitutions <- list(matrix(c(0,0,0), 1, 3)); tips <- 1; print("A total extinction happened!")

								}

							}

							#print(paste("extinction!"))
							
							next

						}
						
						subsevent <- rbinom(seqlen, 1, (mu.new / (1 - Pext)))
						
						if(sum(subsevent) >= 1){
							#print("substitution!")

							# A substitution event occurs. To do this, we bind a row to the substitutions matrix; the columns of this matrix are a site, an initial base, and an end base. If a substitution has already been done in a site, the end base is replaced according to the most recent base.
							
							substsite <- which(subsevent == 1)
							
							for(j in substsite){
							
								if(j %in% substitutions[[i]][, 1]){
								
									substitutions[[i]][which(substitutions[[i]][, 1] == j), 2] <- substitutions[[i]][which(substitutions[[i]][, 1] == j), 3]
								
									substitutions[[i]][which(substitutions[[i]][, 1] == j), 3] <- rownames(subsPmat)[which(rmultinom(1, 1, subsPmat[, substitutions[[i]][which(substitutions[[i]][, 1] == j), 2]]) == 1)]
							
								} else {
							
									substitutions[[i]] <- rbind(substitutions[[i]], c(j, sequence[j], rownames(subsPmat)[which(rmultinom(1, 1, subsPmat[, sequence[j]]) == 1)]))
							
								}
							
							}

							# Finally we add the time period:
							edgetable[i, 3] <- edgetable[i, 3] + dt

							# And the trait or mu and spr must evolve:
							if(direct == F){

								edgetable[i, 4] <- rt + edgetable[i, 4] + rnorm(1, 0, sqrt(2 * D * dt))

							} else {

								propevol <- rnorm(1, 1, 0.05)

								edgetable[i, 4] <- edgetable[i, 4] * propevol

								edgetable[i, 5] <- edgetable[i, 5] * propevol

							}

						}
						
						specevent <- rbinom(1, 1, (spr.new / (1 - Pext)))
						
						if(specevent == 1){
							#print("speciation!")

							# A speciation event occurs and a time step is NOT added:

							if(direct == F){

						   		newbr1 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 1), 0, (rt + edgetable[i, 4] + rnorm(1, 0, sqrt(2 * D * dt))))

						   		newbr2 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 2), 0, (rt + edgetable[i, 4] + rnorm(1, 0, sqrt(2 * D * dt))))

						   	} else {

						   		propevol <- rnorm(1, 1, 0.05)

						   		newbr1 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 1), 0, edgetable[i, 4] * propevol, edgetable[i, 5] * propevol)

						   		newbr2 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 2), 0, edgetable[i, 4] * propevol, edgetable[i, 5] * propevol)

						   	}

						   	edgetable <- rbind(edgetable, newbr1, newbr2)
						   	
						   	substitutions <- append(substitutions, substitutions[length(substitutions)])
						   	
						   	substitutions <- append(substitutions, substitutions[length(substitutions)])

						}
						
						

					} else {

						# Even if no event occurs, we must add a time step:
						edgetable[i, 3] <- edgetable[i, 3] + dt

						# And the trait or mu and spr must evolve:
						if(direct == F){

							edgetable[i, 4] <- rt + edgetable[i, 4] + rnorm(1, 0, sqrt(2 * D * dt))

						} else {

							propevol <- rnorm(1, 1, 0.05)

							edgetable[i, 4] <- edgetable[i, 4] * propevol

							edgetable[i, 5] <- edgetable[i, 5] * propevol

						}

					}

		    	}


		    }

	}


	rownames(edgetable) <- NULL

    colnames(edgetable) <- NULL

    print(time)
    #print(c(extantnodesandtips = length((1:nrow(edgetable))[-extinct]), allnodesandtips = length(1:nrow(edgetable))))

	### The following section takes the object edgetable and spits out two phylogenies, one with branch lengths in terms of time, with any desired total age specified with the argument age, and the other in terms of substitutions.

	simtrtable <- edgetable

    # Chop off the root:
    simtrtable2 <- simtrtable[2:length(simtrtable[,1]), ]
    substitutions2 <- lapply(substitutions, function(x) x <- x[2:nrow(x),])
    substitutions2 <- substitutions2[2:length(substitutions)]
    allsubs <- substitutions2
    
    # Remove all non-extant subsitution data:
    
    for(i in 1:length(substitutions2)){
    	#if(i %in% extinct){
    	#	substitutions2[[i]] <- NA
    	#}
    	#if(simtrtable2[i, 2] %in% simtrtable2[, 1]){
    	#	substitutions2[[i]] <- NA
    	#}
    }

    # Find the number of tips:
    tips <- length(simtrtable2[, 2][which(!simtrtable2[, 2] %in% simtrtable2[, 1])])
    #print(c(alltips = tips, extinctions = length(extinct)))

    # Fix the order so the smallest number is the root:
    fixmat <- rbind(c(0:max(simtrtable2[,2])), c(max(simtrtable2[,2]):0))
    for(i in 1:length(simtrtable2[,1])){
		initval <- simtrtable2[i, 1]
		simtrtable2[i, 1] <- fixmat[2,][which(fixmat[1,] == initval)]
		initval <- simtrtable2[i, 2]
		simtrtable2[i, 2] <- fixmat[2,][which(fixmat[1,] == initval)]
    }

    # Order the branches so the edge object and substitutions are "postorder":
    simtrtabord <- simtrtable2[order(simtrtable2[,1], decreasing = T), ]
    subsord <- substitutions2[order(simtrtable2[,1], decreasing = T)]
    
    # Add one to all values to avoid having a zero:
    simtrtabord[, 1] <- simtrtabord[, 1] + 1

    simtrtabord[, 2] <- simtrtabord[, 2] + 1

    # Fix the values in $edge so that all tips are named 1:n
    currtips <- simtrtabord[,2][which(!simtrtabord[,2] %in% simtrtabord[,1])]

    currtips <- currtips[which(currtips > tips)]

    wrongbrs <- unique(simtrtabord[, 1][which(simtrtabord[, 1] <= tips)])

    for(i in 1:length(currtips)){
		wrongparloc <- which(simtrtabord[, 1] == wrongbrs[i])
		wrongdauloc <- which(simtrtabord[, 2] == currtips[i])
		simtrtabord[, 1][wrongparloc] <- currtips[i]
		simtrtabord[, 2][which(simtrtabord[, 2] == wrongbrs[i])] <- currtips[i]
		simtrtabord[, 2][wrongdauloc] <- wrongbrs[i]
    }

    # Order the branches so the edge and substitutions objects are "postorder" again:
    simtrtabord <- simtrtabord[order(simtrtabord[, 1], decreasing = T), ]
    subsord <- subsord[order(simtrtabord[, 1], decreasing = T)]
    names(subsord) <- simtrtabord[, 2]

    # Extract the edge object and each branch length vector from the table:
    simtredge <- simtrtabord[, 1:2]

    brlentime <- simtrtabord[, 3]

    brlensubs <- simtrtabord[, 4]

    timephylo <- reorder.phylo(rtree(tips), "postorder")

    timephylo$edge <- simtredge

    timephylo$edge.length <- simtrtabord[, 3]
    
    timephylo$tip.label <- 1:tips

    timephylo <- read.tree(text = write.tree(timephylo))

    # Now we modify the age of the phylogeny to be the one given by the argument age.

	if(is.ultrametric(timephylo) == TRUE){

		brlen <- vector()
		for(j in 1:length(timephylo$edge.length)){
			brlen[j] <- (age / max(branching.times(timephylo))) * timephylo$edge.length[j]
		}
		timephylo$edge.length <- brlen

	} else {

		brlen <- vector()
		extantedgelen <- max(timephylo$edge.length[as.vector(which(timephylo$edge[,1] == as.numeric(names(which(branching.times(timephylo) == min(branching.times(timephylo)))))))])
		addedval <- abs(min(branching.times(timephylo))) + extantedgelen
		for(i in 1:length(timephylo$edge.length)){
			brlen[i] <- (age / (max(branching.times(timephylo)) + addedval)) * timephylo$edge.length[i]
		}
		timephylo$edge.length <- brlen

	}

    timephylcut <- drop.fossil(timephylo)
    
    # Now we create the alignment from the substitutions at the tips:
    print(class(subsord));print(dim(subsord[[1]]))
    alignment <- matrix(sequence, nrow = 1, ncol = seqlen, byrow = T)
    
    for(i in 1:length(subsord)){
    	if(is.na(subsord[i]) == F){
    		alignment <- rbind(alignment, sequence)
    		for(j in 1:nrow(subsord[[i]])){
    			alignment[nrow(alignment), as.numeric(subsord[[i]][j, 1])] <- subsord[[i]][j, 3]
    		}
    	}
    }
    print(dim(alignment))
    alignment <- as.DNAbin(alignment[2:nrow(alignment),])
    
    # Finally the function returns a list with the following objects: The complete chronogram, the extant chronogram, and the DNA sequence alignment.

    return(list(timephyloFULL = timephylo, timephylo = timephylcut, alignment = alignment, allsubs = allsubs, subsord = subsord))


}

