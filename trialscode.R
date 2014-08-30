bananas <- list()
pears <- list()
apples <- list()
corban <- vector()
corpea <- vector()
corapp <- vector()
for(i in 1:10){
banana <- try(tr.mu.sp(seqlen = 2000, meanmu = -6, covariance = 25e-6))
try(image(banana[[4]]))
tr <- try(nj(dist.dna(banana[[4]])))
try(print(paste("Distance is", dist.topo(banana[[2]], tr))))
corb <- try(pathnode(tr))
corban[i] <- try(cor(corb[[1]], corb[[2]]))
bananas[[i]] <- try(findroot(banana[[2]], nj(dist.dna(banana[[4]]))))
pear <- try(tr.mu.sp(seqlen = 2000, direct = T, meanmu = -6, covariance = 25e-6))
try(image(pear[[4]]))
tr <- try(nj(dist.dna(pear[[4]])))
try(print(paste("Distance is", dist.topo(pear[[2]], nj(dist.dna(pear[[4]]))))))
corp <- try(pathnode(tr))
corpea[i] <- try(cor(corp[[1]], corp[[2]]))
pears[[i]] <- try(findroot(pear[[2]], nj(dist.dna(pear[[4]]))))
apple <- try(tr.mu.sp(seqlen = 2000, regcoefmu = 0, regcoefspr = 0, meanmu = -6, covariance = 25e-6))
try(image(apple[[4]]))
tr <- try(nj(dist.dna(apple[[4]])))
try(print(paste("Distance is", dist.topo(apple[[2]], nj(dist.dna(apple[[4]]))))))
cora <- try(pathnode(tr))
corapp[i] <- try(cor(cora[[1]], cora[[2]]))
apples[[i]] <- try(findroot(apple[[2]], nj(dist.dna(apple[[4]]))))
}
cortable <- data.frame(corban, corpea, corapp)
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