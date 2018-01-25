# example from: https://msprime.readthedocs.io/en/stable/tutorial.html#replication
# modified to get sequence variants

import msprime
import numpy as np

def segregating_sites_example(n, theta, num_replicates):
    S = np.zeros(num_replicates)
    replicates = msprime.simulate(
                                  sample_size=n,
                                  mutation_rate=theta / 4,
                                  num_replicates=num_replicates)
    for j, tree_sequence in enumerate(replicates):
        S[j] = tree_sequence.get_num_mutations()
        for v in tree_sequence.variants():
            print(v.index, v.position, v.genotypes, sep="\t")

    # Now, calculate the analytical predictions
    S_mean_a = np.sum(1 / np.arange(1, n)) * theta
    S_var_a = (
           theta * np.sum(1 / np.arange(1, n)) +
           theta**2 * np.sum(1 / np.arange(1, n)**2))
    print("              mean              variance")
    print("Observed      {}\t\t{}".format(np.mean(S), np.var(S)))
    print("Analytical    {:.5f}\t\t{:.5f}".format(S_mean_a, S_var_a))

segregating_sites_example(10, 5, 1000)
