require(phangorn)

test.link <- function(phy, prop = 1){
	  nodes <- (length(phy$tip.label)+1:length(phy$tip.label)-1)
	  testnodes <- sample(nodes, length(nodes) * prop)
	  linknodes <- 0

	  for(i in 1:length(testnodes)){
	  	dec1 <- Descendants(phy, Children(phy, testnodes[i])[1])[[1]]
		dec2 <- Descendants(phy, Children(phy, testnodes[i])[2])[[1]]
		div1 <- length(dec1)
		div2 <- length(dec2)
		if(length(dec1) != 1){
			tip1 <- sample(dec1, 1)
			tip2 <- sample(dec2, 1)
		} else {
			tip1 <- dec1
			tip2 <- dec2
		}
		tipname1 <- phy$tip.label[tip1]
		tipname2 <- phy$tip.label[tip2]
		phy <- drop.tip(phy, dec1[-which(dec1 == tip1)])
		phy <- drop.tip(phy, dec1[-which(dec2 == tip2)])
		tp1 <- which(phy$tip.label == tipname1)
		tp2 <- which(phy$tip.label == tipname2)
		br1 <- which(phy$edge[,2] == tp1)
		br2 <- which(phy$edge[,2] == tp2)
		len1 <- phy$edge.length[br1]
		len2 <- phy$edge.length[br2]
		divdif <- div1 - div2
		lendif <- len1 - len2
		if(divdif == 0){ next }
		if(sign(divdif) == sign(lendif)){
				linknodes <- append(linknodes, 1)
		}
	}
	proplink <- sum(linknodes) / (length(nodes) * prop)
}