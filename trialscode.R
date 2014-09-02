bananas <- list()
pears <- list()
apples <- list()
corban <- vector()
corpea <- vector()
corapp <- vector()
for(i in 1:100){
      banana <- try(tr.mu.sp(seqlen = 2000, meanmu = -6, covariance = 25e-6))
      corban[i] <- try(blabla)
      while(class(corban[i]) == "character"){ 
      if(class(banana) == "try-error") banana <- try(tr.mu.sp(seqlen = 2000, meanmu = -6, covariance = 25e-6))
      try(image(banana[[2]]))
      bananas[[i]] <- try(findroot(banana[[1]], nj(dist.dna(banana[[2]]))))
      try(print(paste("Distance is", dist.topo(banana[[1]], bananas[[i]]))))
      corlink <- try(pathnode(bananas[[i]]))
      nodes <- try(as.numeric(corlink[[2]]))
      rtt <- try(as.numeric(corlink[[1]]))
      #model <- nls(nodes ~ d*rtt, start = list(d = 0))
      #lines(rtt, predict(model))
      #corban[i] <- summary(model)$coef[[1]]
      corban[i] <- try(as.numeric(cor(rtt, nodes)))
      corban <- try(as.numeric(corban))
      }
      pear <- try(tr.mu.sp(seqlen = 2000, direct = T, meanmu = -6, covariance = 25e-6))
      corpea[i] <- try(blabla)
      while(class(corpea[i]) == "character"){ 
      if(class(pear) == "try-error") pear <- try(tr.mu.sp(seqlen = 2000, direct = T, meanmu = -6, covariance = 25e-6))
      try(image(pear[[2]]))
      pears[[i]] <- try(findroot(pear[[1]], nj(dist.dna(pear[[2]]))))
      try(print(paste("Distance is", dist.topo(pear[[1]], pears[[i]]))))
      corlink <- try(pathnode(pears[[i]]))
      nodes <- try(as.numeric(corlink[[2]]))
      rtt <- try(as.numeric(corlink[[1]]))
      #model <- nls(nodes ~ d*rtt, start = list(d = 0))
      #lines(rtt, predict(model))
      #corpea[i] <- summary(model)$coef[[1]]
      corpea[i] <- try(as.numeric(cor(rtt, nodes)))
      corpea <- try(as.numeric(corpea))
      }
      apple <- try(tr.mu.sp(seqlen = 2000, regcoefmu = 0, regcoefspr = 0, meanmu = -6, covariance = 25e-6))
      corapp[i] <- try(blabla)
      while(class(corapp[i]) == "character"){ 
      if(class(apple) == "try-error") apple <- try(tr.mu.sp(seqlen = 2000, regcoefmu = 0, regcoefspr = 0, meanmu = -6, covariance = 25e-6))
      try(image(apple[[2]]))
      apples[[i]] <- try(findroot(apple[[1]], nj(dist.dna(apple[[2]]))))
      try(print(paste("Distance is", dist.topo(apple[[1]], apples[[i]]))))
      corlink <- try(pathnode(apples[[i]]))
      nodes <- try(as.numeric(corlink[[2]]))
      rtt <- try(as.numeric(corlink[[1]]))
      #model <- nls(nodes ~ d*rtt, start = list(d = 0))
      #lines(rtt, predict(model))
      #corapp[i] <- summary(model)$coef[[1]]
      corapp[i] <- try(as.numeric(cor(rtt, nodes)))
      corapp <- try(as.numeric(corapp))
      }
      cortable <- data.frame(as.numeric(corban), as.numeric(corpea), as.numeric(corapp))      
      print(cortable)
      boxplot(cortable)
      print(paste("THERE HAVE BEEN", i, "SIMULATIONS THUS FAR"))
}
cortable <- data.frame(as.numeric(corban), as.numeric(corpea), as.numeric(corapp))
print(cortable)
boxplot(cortable)
dev.new()
par(mfrow = c(3, 10), mar = c(1,1,1,1))
for(i in 1:10){
try(plot(ladderize(bananas[[i]]), show.tip.label = F))
}
for(i in 1:10){
try(plot(ladderize(pears[[i]]), show.tip.label = F))
}
for(i in 1:10){
try(plot(ladderize(apples[[i]]), show.tip.label = F))
}