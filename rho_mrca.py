# example rho calculation
# important: rho implemented for MRCA only

import msprime
import numpy as np

import matplotlib.pyplot as plot

# model parameters
random_seed        = 10
num_replicates     = 1000
sample_size        = 20
Ne                 = 1e4
length             = 5e3
mutation_rate      = 2e-8
recombination_rate = 0

replicates = msprime.simulate(
    sample_size=sample_size, 
    Ne=Ne, 
    length=length, 
    mutation_rate=mutation_rate, 
    recombination_rate=recombination_rate,
    random_seed=random_seed,
    num_replicates=num_replicates)

def rho_mrca(replicates):
		
	means = []
	
	for j, tree_sequence in enumerate(replicates):
		
		# transpose genotype matrix
		variants = tree_sequence.genotype_matrix().transpose()
		
		# generate array of lineage elements
		l = variants.sum(axis=1)
		
		# print mean, variance, standard deviation, lower 95% CI, upper 95% CI
		# using ddof N-1 correction (matches R default)
		print(j+1, np.mean(l), np.var(l, ddof = 1), np.std(l, ddof = 1), np.percentile(l, 2.5), np.percentile(l, 97.5), sep='\t')
		
		# save mean value to array
		means.append(np.mean(l))
	
	return(means)
		
means = rho_mrca(replicates)


# plot rho
plot.hist(means,density=1, bins=500) 
plot.xlabel('Rho')
plot.ylabel('Density')

plot.show()
