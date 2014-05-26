require(phangorn)

test.link <- function(phy, prop = 1){
	  nodes <- (length(phy$tip.label)+1):(length(phy$tip.label)*2-1)
	  testnodes <- sample(nodes, length(nodes) * prop)
	  linknodes <- 0

	  for(i in 1:length(testnodes)){
	  	dec1 <- Descendants(phy, Children(phy, testnodes[i])[1])[[1]]
		dec2 <- Descendants(phy, Children(phy, testnodes[i])[2])[[1]]
		div1 <- length(dec1)
		div2 <- length(dec2)
		if(div1 > 1){ tip1 <- sample(dec1, 1) } else {
			tip1 <- dec1
		}
                if(div2 > 1){ tip2 <- sample(dec2, 1) } else {
                        tip2 <- dec2
                }
		tipname1 <- phy$tip.label[tip1]
		tipname2 <- phy$tip.label[tip2]
                if(div1 > 1){
                    phy1 <- drop.tip(phy, dec1[-which(dec1 == tip1)])
                }
                if(div2 > 1){
                    phy1 <- drop.tip(phy, dec2[-which(dec2 == tip2)])
                }
		tp1 <- which(phy$tip.label == tipname1)
		tp2 <- which(phy$tip.label == tipname2)
		br1 <- which(phy$edge[,2] == tp1)
		br2 <- which(phy$edge[,2] == tp2)
		len1 <- phy$edge.length[br1]
		len2 <- phy$edge.length[br2]
		divdif <- div1 - div2
		lendif <- len1 - len2
                #print(paste(i, divdif, lendif))

                if(divdif == 0 || lendif == 0){ next }
		if(sign(divdif) == sign(lendif)){
				linknodes <- append(linknodes, 1)
		}
	}
	proplink <- sum(linknodes) / length(testnodes)
        return(proplink)
}

### The following is code to test the quality of simulations using tr.mu.sp

nolinklink <- 0
trlinklink <- 0
dirlinklink <- 0
nolinkgamma <- 0
trlinkgamma <- 0
dirlinkgamma <- 0
for(i in 1:100){
trnolink <- simtr2seq(tr.mu.sp(regcoefmu = 0, regcoefspr = 0, stepsize = 0.001))
trtrlink <- simtr2seq(tr.mu.sp(stepsize = 0.001))
trdirlink <- simtr2seq(tr.mu.sp(molerror = 0, sprerror = 0, stepsize = 0))
nolinklink <- append(nolinklink, test.link(trnolink[[4]]))
trlinklink <- append(trlinklink, test.link(trtrlink[[4]]))
dirlinklink <- append(dirlinklink, test.link(trdirlink[[4]]))
nolinkgamma <- append(nolinkgamma, gammaStat(trnolink[[3]]))
trlinkgamma <- append(trlinkgamma, gammaStat(trtrlink[[3]]))
dirlinkgamma <- append(dirlinkgamma, gammaStat(trdirlink[[3]]))
}
