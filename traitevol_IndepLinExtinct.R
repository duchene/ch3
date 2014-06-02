
require(phangorn)
require(MASS)

# The following function simulates phylogenetic structures (not yet of class phylo) according to the relationsips between a trait and diversification rate and the same trait and the probability of a substitution occurring.

# Arguments:	stepsize is the size of each time step in time (it influences the probability of an event of speciation, exinction, or substitution occurring).	branchstop is number of branches desired and stops the simulation. 	seqlen is the length of the genetic sequence to be generated.	traitstart is the initial value of the trait.			trait.r is the rate of change of the trait, which is adjusted when a trend is desired. 		FUNspr and FUNmu are the functions that define the relationship between the trait and the probability of bifurcation and substitution respectively.	Pext is the constant background probability of extinction.	D is the variance of trait evolution, which is taken as being constant.	molerror and sprerror are the error to be introduced to each of speciation and mutation probability at each step.


tr.mu.sp <- function(stepsize = 0.01, branchstop = 200, seqlen = 2000, traitstart = 50, trait.r = 0, regcoefmu = 0.01, regcoefspr = 0.02, Pext = 0.01, Dsd = 0.001, molerror = 0.001, sprerror = 0.01, direct = F, covariance = 0.0000099, meanmu = 0.009, meanspr = 0.2){

# The following is a matrix where the columns are: the parent node, the daughter node, the branch length in time, the branch length in substitutions, and etiher the trait value or the means for mu and spr.
	
	if(direct == F){
	
		edgetable <- matrix(data = c(0, 1, 0, 0, traitstart), nrow = 1, ncol = 5)
		
	} else {
		
		edgetable <- matrix(data = c(0, 1, 0, 0, meanmu, meanspr), nrow = 1, ncol = 6)
		
	}
	 
	extinct <- 1
	 
	time <- 1 
	 
	b1 <- regcoefmu
		    		
	b2 <- regcoefspr
	
	muvar <- molerror^2
	
	spvar <- sprerror^2
	 
	 while(length((1:nrow(edgetable))[if(length(extinct) > 1){ -extinct } else { 1:nrow(edgetable) }]) < branchstop){

	 	    # Here we are setting constant time step size. Otherwise, use rexp().
	 	   
		    dt <- stepsize
        	rt <- trait.r * dt
        	
        	# The following restarts the simulation if the age of the tree exceeding 100.
        	time <- time + dt
        	if(time > 501){ time <- 1; extinct <- 1; 
        		
        		if(direct == F){
        		
        			edgetable <- matrix(data = c(0, 1, 0, 0, traitstart), nrow = 1, ncol = 5); print("A phylogeny grew too old!")
        			
        		} else {
        			
        			edgetable <- matrix(data = c(0, 1, 0, 0, meanmu, meanspr), nrow = 1, ncol = 6); print("A phylogeny grew too old")
        			
        		}
        		
        	}
        	

		    # If interested on the rate of evolution of the trait, we should set D differently, perhaps with dexp().

		    for(i in (1:nrow(edgetable))[if(length(extinct) > 1){ -extinct } else { 1:nrow(edgetable) }]){
		    
		    	D <- rexp(1, Dsd)

		    	if(!(edgetable[i, 2] %in% edgetable[, 1])){ # Grab a tip lineage
					
					e1 <- rnorm(1, 0, molerror)
		    		
					e2 <- rnorm(1, 0, sprerror)

					FUNspr <- function(x) abs(0.3 * (1 - exp(-b2 * x)) + 0.02 + e2)

					FUNmu <- function(x) abs(0.02 * (1 - exp(-b1 * x)) + 0.001 + e1)
	
					if(b2 == 0){
						FUNspr <- function(x) abs(0.3 * (1 - exp(-b2 * x)) + meanspr + e2)
					}
					if(b1 == 0){
						FUNmu <- function(x) abs(0.02 * (1 - exp(-b1 * x)) + meanmu + e1)
					}
					
					if(direct == F){
					
						# The following makes sure that trait values do not become lower than roughly 2.

						if(edgetable[i, 5] < 2) edgetable[i, 5] <- 2 + rnorm(1, 0, sqrt(2 * D * dt))

		    			spr.new <- FUNspr(edgetable[i, 5])
		    		
		    			mu.new <- FUNmu(edgetable[i, 5]) * seqlen
		    		
		    		} else {
		    			
		    			if(edgetable[i, 5] < 0.001) edgetable[i, 5] <- 0.001 + rnorm(1, 0, 0.0001)
		    			
		    			if(edgetable[i, 6] < 0.01) edgetable[i, 6] <- 0.01 + rnorm(1, 0, 0.001)
		    			
		    			rand <- mvrnorm(1, c(edgetable[i, 5], edgetable[i, 6]), matrix(c(muvar, covariance, covariance, spvar), 2, 2))
		    			
		    			while(any(rand < 0)) rand <- mvrnorm(1, c(edgetable[i, 5], edgetable[i, 6]), matrix(c(muvar, covariance, covariance, spvar), 2, 2))
		    			
		    			spr.new <- rand[2]
		    			
		    			mu.new <- rand[1] * seqlen
		    		
		    		}
		    		
		    		# The following is the total probability of an event occurring, which includes twice the equation for the rate of an event occurring.

		    		totalP <- ((spr.new + mu.new + Pext) * dt) * exp(-((spr.new + mu.new + Pext) * dt))

						# The following condition determines whether and event occurs.

		    		if(runif(1) <= totalP){

		    			specP <- spr.new / (spr.new + mu.new + Pext)
		    				
		    			extP <- Pext / (spr.new + mu.new + Pext)
		    				
		    			muP <- mu.new / (spr.new + mu.new + Pext)

		    			# The following segment determines which event occurs.
		    				
		    			event <- which(rmultinom(1, 1, c(specP, extP, muP)) == 1)

						if(event == 1){
							#print("speciation!")

							# A speciation event occurs and a time step is added:
							
							if(direct == F){
							
						   		newbr1 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 1), dt, 0, (rt + edgetable[i, 5] + rnorm(1, 0, sqrt(2 * D * dt))))

						   		newbr2 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 2), dt, 0, (rt + edgetable[i, 5] + rnorm(1, 0, sqrt(2 * D * dt))))
						   	
						   	} else {
						   		
						   		propevol <- rnorm(1, 1, 0.05)
						   	
						   		newbr1 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 1), dt, 0, edgetable[i, 5] * propevol, edgetable[i, 6] * propevol)

						   		newbr2 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 2), dt, 0, edgetable[i, 5] * propevol, edgetable[i, 6] * propevol)
						   	
						   	}

						   	edgetable <- rbind(edgetable, newbr1, newbr2)

						} else if(event == 2){
							
							# An extinction event occurs and a time step is added:
							
							extinct <- append(extinct, i)
							
							if(length(extinct) > length(!(edgetable[, 2] %in% edgetable[, 1]))) { 
							
								if(direct == F){
							
									edgetable <- matrix(data = c(0, 1, 0, 0, traitstart), nrow = 1, ncol = 5); extinct <- 1 ; time <- 1; print("A total extinction happened!")
							
								} else {
							
									edgetable <- matrix(data = c(0, 1, 0, 0, meanmu, meanspr), nrow = 1, ncol = 6); extinct <- 1 ; time <- 1; print("A total extinction happened!")
							
								}
							
							}
							
							#print(paste("extinction!"))
							
						} else if(event == 3){
							#print("mutation!")

							# A mutation event occurs:
							edgetable[i, 4] <- edgetable[i, 4] + 1

							# Finally we add the time period:
							edgetable[i, 3] <- edgetable[i, 3] + dt

							# And the trait or mu and spr must evolve:
							if(direct == F){
								
								edgetable[i, 5] <- rt + edgetable[i, 5] + rnorm(1, 0, sqrt(2 * D * dt))
								
							} else {
							
								propevol <- rnorm(1, 1, 0.05)
								
								edgetable[i, 5] <- edgetable[i, 5] * propevol
								
								edgetable[i, 6] <- edgetable[i, 6] * propevol
							
							}

						}

					} else {

						# Even if no event occurs, we must add a time step:
						edgetable[i, 3] <- edgetable[i, 3] + dt

						# And the trait or mu and spr must evolve:
						if(direct == F){
						
							edgetable[i, 5] <- rt + edgetable[i, 5] + rnorm(1, 0, sqrt(2 * D * dt))
						
						} else {
						
							propevol <- rnorm(1, 1, 0.05)
							
							edgetable[i, 5] <- edgetable[i, 5] * propevol
								
							edgetable[i, 6] <- edgetable[i, 6] * propevol
						
						}

					}

		    	}


		    }

	}


	rownames(edgetable) <- NULL

    colnames(edgetable) <- NULL
    
    edgetable[, 4] <- sapply(edgetable[, 4], function(x) if(x == 0){ x <- 0.0001 } else { x <- x })
    
    print(time)

	return(edgetable)	

}


### The following function takes the object created by tr.mu.sp and spits out two phylogenies, one with branch lengths in terms of time, with any desired total age specified with the argument age, and the other in terms of substitutions.


simtr2seq <- function(simtrtable, q = rep(1, 6), seqlen = 2000, age = 50){

    # Chop off the root:
    simtrtable2 <- simtrtable[2:length(simtrtable[,1]), ]

    # Find the number of tips:
    tips <- length(simtrtable2[, 2][which(!simtrtable2[, 2] %in% simtrtable2[, 1])])

    # Fix the order so the smallest number is the root:
    fixmat <- rbind(c(0:max(simtrtable2[,2])), c(max(simtrtable2[,2]):0))
    for(i in 1:length(simtrtable2[,1])){
	initval <- simtrtable2[i, 1]
	simtrtable2[i, 1] <- fixmat[2,][which(fixmat[1,] == initval)]
	initval <- simtrtable2[i, 2]
	simtrtable2[i, 2] <- fixmat[2,][which(fixmat[1,] == initval)]
    }

    # Order the branches so the edge object is "postorder":
    simtrtabord <- simtrtable2[order(simtrtable2[,1], decreasing = T), ]
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

    # Order the branches so the edge object is "postorder" again:
    simtrtabord <- simtrtabord[order(simtrtabord[, 1], decreasing = T), ]

    # Extract the edge object and each branch length vector from the table:
    simtredge <- simtrtabord[, 1:2]

    brlentime <- simtrtabord[, 3]

    brlensubs <- simtrtabord[, 4]

    timephylo <- reorder.phylo(rtree(tips), "postorder")

    timephylo$edge <- simtredge

    timephylo$edge.length <- simtrtabord[, 3]

    subsphylo <- timephylo

    subsphylo$edge.length <- simtrtabord[, 4] / (10 * max(simtrtabord[, 4]))
    
    subsphylo$edge.length[which(subsphylo$edge == 0)] <- 1 / seqlen

    timephylo <- read.tree(text = write.tree(timephylo))

    subsphylo <- read.tree(text = write.tree(subsphylo))
    
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
    
    # Finally the function returns a list with the following objects: The complete chronogram, the complete phylogram, the extant chronogram, the extant phylogram, and the DNA sequence alignment.
    
    timephylcut <- drop.fossil(timephylo)
    
    subsphylcut <- drop.tip(subsphylo, setdiff(subsphylo$tip.label, timephylcut$tip.label))
    
    DNAalignment <- as.DNAbin(simSeq(subsphylcut, l = seqlen, Q = q, ancestral = F))

    return(list(timephyloFULL = timephylo, subsphyloFULL = subsphylo, timephylo = timephylcut, subsphylo = subsphylcut, alignment = DNAalignment))
    
    
}




    
