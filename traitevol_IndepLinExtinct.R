
require(phangorn)
require(MASS)
require(Matrix)

# The following function simulates phylogenetic structures (not yet of class phylo) according to the relationsips between a trait and diversification rate and the same trait and the probability of a substitution occurring.

# Arguments:	stepsize is the size of each time step in time (it influences the probability of an event of speciation, exinction, or substitution occurring).	branchstop is number of branches desired and stops the simulation. 	seqlen is the length of the genetic sequence to be generated.	traitstart is the initial value of the trait.			trait.r is the rate of change of the trait, which is adjusted when a trend is desired. 		FUNspr and FUNmu are the functions that define the relationship between the trait and the probability of bifurcation and substitution respectively.	Pext is the constant background probability of extinction.	D is the variance of trait evolution, which is taken as being constant.	molerror and sprerror are the error to be introduced to each of speciation and mutation probability at each step.


tr.mu.sp <- function(stepsize = 0.1, branchstop = 200, seqlen = 2000, traitstart = 50, trait.r = 0, regcoefmu = 0.58, regcoefspr = 1.15, Pext = 0.01, D = 1, molerror = 0.0005, sprerror = 0.05, direct = F, covariance = 25e-6, meanmu = -4.6, meanspr = -2.3, q = 0, age = 0){

# The following is a matrix where the columns are: the parent node, the daughter node, the branch length in time, and etiher the trait value or the means for mu and spr.

	if(direct == F){

		edgetable <- matrix(data = c(0, 1, 0, traitstart), nrow = 1, ncol = 4)

	} else {

		edgetable <- matrix(data = c(0, 1, 0, meanmu, meanspr), nrow = 1, ncol = 5)

	}

	if(q == 0){
	     q <- matrix(c(-7.0921, 1.3472, 4.8145, 0.9304, 1.3472, -8.155, 1.2491, 5.5587, 4.8145, 1.2491, -7.0636, 1, 0.9304, 5.5587, 1, -7.4891), 4, 4, byrow = T) # From placental mammals. In Murphy et al. 2001 in Science
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

	while(((tips - (length(extinct) - 1)) * 2) < branchstop){

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
        	
		    # If interested on the rate of evolution of the trait, we should set D differently, perhaps with dexp().

		    for(i in (1:nrow(edgetable))[if(length(extinct) > 1){ -extinct } else { 1:nrow(edgetable) }]){

		    	if(!(edgetable[i, 2] %in% edgetable[, 1])){ # Grab a tip lineage

					e1 <- rnorm(1, 0, molerror)

					e2 <- rnorm(1, 0, sprerror)

					FUNspr <- function(x) b2 * log(x) - 6.9 + e2

					FUNmu <- function(x) b1 * log(x) - 6.9 + e1

					if(b2 == 0){
						FUNspr <- function(x) meanspr #+ (e2 / 10)
					}
					if(b1 == 0){
						FUNmu <- function(x) meanmu #+ (e1 / 10)
					}

					if(direct == F){

						# The following makes sure that trait values do not become lower than roughly 2.

						if(edgetable[i, 4] < 2) edgetable[i, 4] <- 2 + rnorm(1, 0, (D * dt))
						if(edgetable[i, 4] > 100) edgetable[i, 4] <- 100 + rnorm(1, 0, (D * dt))

		    			spr.new0 <- suppressWarnings(FUNspr(edgetable[i, 4]))

		    			mu.new0 <- suppressWarnings(FUNmu(edgetable[i, 4]))
		    			
		    			spr.new <- exp(spr.new0)
		    			
		    			mu.new <- exp(mu.new0)
		    			
		    			#print(edgetable[i, 4])
		    			
		    			#print(c(spr.new0, spr.new, mu.new0, mu.new))


		    		} else {

		    			rand <- mvrnorm(1, c(edgetable[i, 4], edgetable[i, 5]), matrix(c(muvar, covariance, covariance, spvar), 2, 2))

		    			if(rand[1] > -4.19) rand[1] <- -4.19
		    			
		    			if(rand[2] > -1.6) rand[2] <- -1.6

		    			spr.new <- exp(rand[2])

		    			mu.new <- exp(rand[1])
		    			
		    			#print(c(spr.new, mu.new))

		    		}
		    		
		    		event <- 0
		    			
		    		### The following is the total probability of each event occurring.

	    			totalP <- ((spr.new + mu.new + Pext + (spr.new * mu.new)) * dt) * exp(-((spr.new + mu.new + Pext + (spr.new * mu.new)) * dt))
	    			
	    			#print(totalP)
	    			
	    			# The following determines whether an event occurs.
	    			
	    			if(runif(1) <= totalP){
	    			
	    				specP <- spr.new / (spr.new + (mu.new * seqlen) + Pext + (spr.new * (mu.new * seqlen)))
		    				
		    			extP <- Pext / (spr.new + (mu.new * seqlen) + Pext + (spr.new * (mu.new * seqlen)))
		    				
		    			muP <- (mu.new * seqlen) / (spr.new + (mu.new * seqlen) + Pext + (spr.new * (mu.new * seqlen)))
		    			
		    			muspP <- (spr.new * (mu.new * seqlen)) / (spr.new + (mu.new * seqlen) + Pext + (spr.new * (mu.new * seqlen)))

		    			# The following segment determines which event occurs.
		    				
		    			event <- which(rmultinom(1, 1, c(extP, muP, specP, muspP)) == 1)
		    			
		    		} else {
		    			
		    			event <- 0
		    		
		    		}

					if(event == 1){

						# An extinction event occurs:

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
						
					if(event == 2 | event == 4){
						
					subsevent <- rbinom(seqlen, 1, mu.new)
						
					if(sum(subsevent) >= 1){
						
						#print(paste("substitution!", sum(subsevent)))

						# A substitution event occurs. To do this, we bind a row to the substitutions matrix; the columns of this matrix are a site, an initial base, and an end base. If a substitution has already been done in a site, the end base is replaced according to the most recent base.
							
						for(j in 1:length(subsevent)){
							
							if(subsevent[j] == 1){
							
								if(j %in% substitutions[[i]][, 1]){
								
									substitutions[[i]][which(substitutions[[i]][, 1] == j), 2] <- substitutions[[i]][which(substitutions[[i]][, 1] == j), 3]
								
									substitutions[[i]][which(substitutions[[i]][, 1] == j), 3] <- rownames(subsPmat)[which(rmultinom(1, 1, subsPmat[, substitutions[[i]][which(substitutions[[i]][, 1] == j), 2]]) == 1)]
							
								} else {
							
									substitutions[[i]] <- rbind(substitutions[[i]], c(j, sequence[j], rownames(subsPmat)[which(rmultinom(1, 1, subsPmat[, sequence[j]]) == 1)]))
							
									}
								
								}
							
							}


						}
						
						}
						
						if(event == 3 | event == 4){
							#print("speciation!")

							# A speciation event occurs and a time step is added:

							if(direct == F){

						   		newbr1 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 1), dt, (rt + edgetable[i, 4] + rnorm(1, 0, (D * dt))))

						   		newbr2 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 2), dt, (rt + edgetable[i, 4] + rnorm(1, 0, (D * dt))))

						   	} else {
							        
								# From here onwards, we need both speciation and substitution to evolve together, so we create a propevol value that allows them to evolve by the same proportion every time.

								propevol <- rnorm(1, 0, 0.05)

						   		newbr1 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 1), dt, edgetable[i, 4] + (meanmu * propevol), edgetable[i, 5] + (meanspr * propevol))

						   		newbr2 <- c(edgetable[i, 2], (max(edgetable[, 2]) + 2), dt, edgetable[i, 4] + (meanmu * propevol), edgetable[i, 5] + (meanspr * propevol))

						   	}

						   	edgetable <- rbind(edgetable, newbr1, newbr2)
						   	
						   	substitutions <- append(substitutions, substitutions[i])
						   	
						   	substitutions <- append(substitutions, substitutions[i])
						   	#print(dim(substitutions[[i]]))

						} 
						
						
						# If no speciation occurred, but a substitution did, we still add a time period:
						
						if(event == 2){
							
							edgetable[i, 3] <- edgetable[i, 3] + dt

							# And the trait or mu and spr must evolve:
							if(direct == F){

								edgetable[i, 4] <- rt + edgetable[i, 4] + rnorm(1, 0, (D * dt))

							} else {

								propevol <- rnorm(1, 0, 0.05)

								edgetable[i, 4] <- edgetable[i, 4] + (meanmu * propevol)

								edgetable[i, 5] <- edgetable[i, 5] + (meanspr * propevol)

							}
						}
						

					if(event == 0){

						# Even if no event occurs, we must add a time step:
						edgetable[i, 3] <- edgetable[i, 3] + dt

						# And the trait or mu and spr must evolve:
						if(direct == F){

							edgetable[i, 4] <- rt + edgetable[i, 4] + rnorm(1, 0, (D * dt))

						} else {

							propevol <- rnorm(1, 0, 0.05)

							edgetable[i, 4] <- edgetable[i, 4] + (meanmu * propevol)

							edgetable[i, 5] <- edgetable[i, 5] + (meanspr * propevol)

						}

					}

		    	}
		    	
        	tips <- length(edgetable[, 2][which(!edgetable[, 2] %in% edgetable[, 1])])

		    }

	}

	rownames(edgetable) <- NULL

    colnames(edgetable) <- NULL

    print(time)
    #print(c(tips = tips, extinct = length(extinct) - 1))
    
	### The following section takes the object edgetable and spits out two phylogenies, one with branch lengths in terms of time, with any desired total age specified with the argument age, and the other in terms of substitutions.

	simtrtable <- edgetable

    # Chop off the root:
    simtrtable2 <- simtrtable[2:nrow(simtrtable), ]
    for(i in 1:length(substitutions)){
    	#print(dim(substitutions[[i]]))
    	substitutions[[i]] <- substitutions[[i]][2:nrow(substitutions[[i]]),]
    }
    
    substitutions2 <- substitutions[2:length(substitutions)]
    
    #print(simtrtable2[1:5, ])
    
    # Create the sequences for all nodes and tips
    
    alignment <- matrix(sequence, nrow = 1, ncol = seqlen, byrow = T)
    
    for(i in 1:length(substitutions2)){
    	alignment <- rbind(alignment, sequence)
    	for(j in 1:nrow(substitutions2[[i]])){
    			alignment[nrow(alignment), as.numeric(substitutions2[[i]][j, 1])] <- substitutions2[[i]][j, 3]
    	}
    }
    #print(dim(alignment))
    simtrtable2 <- suppressWarnings(cbind(as.data.frame(simtrtable2), as.data.frame(alignment[2:nrow(alignment), ])))
    #print(class(simtrtable2))
    #print(dim(simtrtable2))
    #print(simtrtable2[1:5, 1:10])
    
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

    # Order the branches so the edge object and substitutions are "postorder":
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

    # Order the branches so the edge and substitutions objects are "postorder" again:
    
    simtrtabord <- simtrtabord[order(simtrtabord[, 1], decreasing = T), ]
    
    # Extract the edge object and each branch length vector from the table:
    simtredge <- simtrtabord[, 1:2]

    brlentime <- simtrtabord[, 3]

    timephylo <- reorder.phylo(rtree(tips), "postorder")

    timephylo$edge <- simtredge

    timephylo$edge.length <- simtrtabord[, 3]
    
    timephylo$tip.label <- 1:tips

    timephylo <- read.tree(text = write.tree(timephylo))

    # Now we modify the age of the phylogeny to be the one given by the argument age.

    if(age > 0){

	if(is.ultrametric(timephylo) == TRUE){

		brlen <- vector()
		for(j in 1:length(timephylo$edge.length)){
			brlen[j] <- (age / max(branching.times(timephylo))) * timephylo$edge.length[j]
		}
		timephylo$edge.length <- brlen

	} else {

		brlen <- vector()
		
		minbrage <- min(branching.times(timephylo))
		
		brtimes <- branching.times(timephylo)
		
		youngbrs <- as.numeric(names(brtimes[which(brtimes == minbrage)]))
		
		#print(youngbrs)
		
		extantedgelen <- 0
		for(i in 1:length(youngbrs)){
		tempextedgelen <- max(timephylo$edge.length[as.vector(which(timephylo$edge[,1] == youngbrs[i] ))])
		extantedgelen <- max(extantedgelen, tempextedgelen)
		}
		addedval <- abs(min(branching.times(timephylo))) + extantedgelen
		for(i in 1:length(timephylo$edge.length)){
			brlen[i] <- (age / (max(branching.times(timephylo)) + addedval)) * timephylo$edge.length[i]
		}
		timephylo$edge.length <- brlen

	}

	}

    timephylcut <- drop.fossil(timephylo)
    
    # Extract alignments of all nodes, all tips, and extant tips.
    
    if(direct == F){
    
    	alignallnodes <-  as.DNAbin(as.matrix(simtrtabord[, 5:ncol(simtrtabord)]))
    	
    	rownames(alignallnodes) <- simtrtabord[, 2]
    
    	alignalltips <- as.DNAbin(as.matrix(simtrtabord[, 5:ncol(simtrtabord)][which(simtrtabord[, 2] %in% 1:tips), ]))
    	
    	rownames(alignalltips) <- simtrtabord[, 2][which(simtrtabord[, 2] %in% 1:tips)]
    	
    	alignextant <- alignalltips[timephylcut$tip.label,]
    
    } else {
    
    	alignallnodes <-  as.DNAbin(as.matrix(simtrtabord[, 6:ncol(simtrtabord)]))
    	
    	rownames(alignallnodes) <- simtrtabord[, 2]
    
    	alignalltips <- as.DNAbin(as.matrix(simtrtabord[, 6:ncol(simtrtabord)][which(simtrtabord[, 2] %in% 1:tips), ]))
    	
    	rownames(alignalltips) <- simtrtabord[, 2][which(simtrtabord[, 2] %in% 1:tips)]
    	
    	alignextant <- alignalltips[timephylcut$tip.label,]
    
    }
    
    # Finally the function returns a list with the following objects: The complete chronogram, the extant chronogram, and the DNA sequence alignment.

    return(list(timephyloFULL = timephylo, timephylo = timephylcut, alltipsDNA = alignalltips, extantDNA = alignextant))


}

