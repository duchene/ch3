for(i in 1:10){
banana <- try(tr.mu.sp(seqlen = 500))
bananas[[i]] <- try(findroot(banana[[2]], nj(dist.dna(banana[[4]]))))
pear <- try(tr.mu.sp(seqlen = 500, direct = T))
pears[[i]] <- try(findroot(pear[[2]], nj(dist.dna(pear[[4]]))))
apple <- try(tr.mu.sp(seqlen = 500, regcoefmu = 0, regcoefspr = 0))
apples[[i]] <- try(findroot(apple[[2]], nj(dist.dna(apple[[4]]))))
}

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