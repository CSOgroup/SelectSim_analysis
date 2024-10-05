# testing null model efficiency

library(tictoc)

# BeRewire version
library(BiRewire)

data(BRCA_binary_matrix)
dim(BRCA_binary_matrix)
# [1]  757 9757
# sample by gene format
brca = t(BRCA_binary_matrix)

gams = list()
gams[[1]] = brca

unit = brca[sample(1:nrow(brca), size = 500), sample(1:ncol(brca), size = 500)]
gams[[2]] = unit

for(i in 1:19){
	temp = cbind(gams[[i+1]],unit)
	temp = matrix(sample(temp), nrow = nrow(unit))
	gams[[i+2]] = temp
}

for(i in 1:length(gams)){
	cat(dim(gams[[i]]),'\n')
}

# gams is a vector of 20 GAMs

# BiRewire execution: compute 100 random permutation for each GAM
br.times = c()
for(i in 1:length(gams)){
	
	max=100*sum(gams[[i]])
	
	tic('test BR:')
	for(k in 1:100){
		rw.brca = birewire.rewire.bipartite(gams[[i]], max.iter = max, verbose=FALSE)
	}
	v = toc()
	
	value = v$toc - v$tic	
	br.times = c(br.times, value)

}


# simulation execution: compute 100 random simulations for each GAM
ss.times = c()
for(i in 1:length(gams)){

	tic('test SS:')
	
	tmb = as.matrix( colSums(gams[[i]])/sum(colSums(gams[[i]])), ncol = 1 )
	gen = as.matrix( rowSums(gams[[i]]), ncol = 1 )

	template = gen %*% t(tmb)
	nvalues = nrow(template)*ncol(template)
	
	for(k in 1:100){
		r = matrix( runif(nvalues, min = 0, max = 1), nrow = nrow(template), ncol = ncol(template))
		S = 1*((template - r) > 0)
	}
	v = toc()
	
	value = v$toc - v$tic	
	ss.times = c(ss.times, value)
	
}
