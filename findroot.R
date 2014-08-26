findroot <- function(sim.tree, est.tree){
	 rootnode <- length(sim.tree$tip.label) + 1
	 outgroupnode <- sample(Children(sim.tree, rootnode), 1)
	 outgroup <- sim.tree$tip.label[Descendants(sim.tree, outgroupnode, type = "tips")[[1]]]
	 newtree <- root(est.tree, outgroup)
	 return(newtree)
}
