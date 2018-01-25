# example Asian-Papuan demography
# originally defined in Sanderson et al. 2015: https://doi.org/10.1534/genetics.115.176842

import msprime
import libsequence.polytable as pt
import libsequence.summstats as sstats
import libsequence.windows as windows

# basic variables
generation_time = 25
length = 1e4
recombination_rate = 2e-8
mutation_rate = 2e-8

# sample sizes
S_asian   = 100
S_papuan  = 100
S_admixed = 100

# effective sizes
N_ancestral = 10500
N_asian     = 2050
N_papuan    = 800
N_admixed   = 1425

# demography times in years
T_ancestral = 50000 / generation_time
T_admixture = 4000 / generation_time

# asian admixture proportion
P_admixture = 0.5

# migration rates
m_AS_PA = 0
m_AS_AD = 0
m_PA_AD = 0

# population IDs: asian=0, papuan=1, admixed=2
population_configurations = [
    msprime.PopulationConfiguration(sample_size=S_asian,   initial_size=N_asian),
    msprime.PopulationConfiguration(sample_size=S_papuan,  initial_size=N_papuan),
    msprime.PopulationConfiguration(sample_size=S_admixed, initial_size=N_admixed)
]

# migration
migration_matrix = [
    [      0, m_AS_PA, m_AS_AD],
    [m_AS_PA,       0, m_PA_AD],
    [m_AS_AD, m_PA_AD,       0],
]

# demography
demographic_events = [
    # Admixed merges to Asian and Papuan
    msprime.MassMigration(time=T_admixture, source=2, destination=0, proportion=P_admixture),
    msprime.MassMigration(time=T_admixture, source=2, destination=1, proportion=1.0),
    msprime.MigrationRateChange(time=T_admixture, rate=0),
    msprime.MigrationRateChange(time=T_admixture, rate=m_AS_PA, matrix_index=(0, 1)),
    msprime.MigrationRateChange(time=T_admixture, rate=m_AS_PA, matrix_index=(1, 0)),
    msprime.PopulationParametersChange(time=T_admixture, initial_size=N_asian, growth_rate=0, population_id=0),
    msprime.PopulationParametersChange(time=T_admixture, initial_size=N_papuan, growth_rate=0, population_id=1),
    #msprime.PopulationParametersChange(time=T_admixture, initial_size=0, growth_rate=0, population_id=2),
    # Asian and Papuan merge
    msprime.MassMigration(time=T_ancestral, source=1, destination=0, proportion=1.0),
    msprime.MigrationRateChange(time=T_ancestral, rate=0),
    msprime.PopulationParametersChange(time=T_ancestral, initial_size=N_ancestral, growth_rate=0, population_id=0),
    #msprime.PopulationParametersChange(time=T_ancestral, initial_size=0, growth_rate=0, population_id=1),
]

# check demography (human readable only)
dp = msprime.DemographyDebugger(
    Ne=N_ancestral,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,
    demographic_events=demographic_events)
dp.print_history()

# simulate demography
sims = msprime.simulate(
    length=length,
    recombination_rate=recombination_rate,
    mutation_rate=mutation_rate,
    Ne=N_ancestral,
    population_configurations=population_configurations,
    migration_matrix=migration_matrix,
    demographic_events=demographic_events)

# get population sample lists
asians  = sims.get_samples(population_id=0)
papuans = sims.get_samples(population_id=1)
admixed = sims.get_samples(population_id=2)
pops = {'Asians': asians, 'Papuans': papuans, 'Admixed': admixed}

# temp: sanity check variants (and population subsetting)
for v in sims.variants():
    print(v.index, v.position, v.genotypes, sep="\t")
print("\n")

for i in pops.keys():
    for v in sims.variants():
        print(v.index, v.position, v.genotypes[pops[i]], sep="\t")
    print("\n")

# use subsetting to access populations
# everyone
sd = pt.SimData([(v.position, v.genotypes) for v in sims.variants(as_bytes=True)])
ps = sstats.PolySIM(sd)
print("Everyone:", ps.tajimasd(),' ',ps.thetaw(),' ',ps.thetapi())

for i in pops.keys():
    sd = pt.SimData([(v.position, v.genotypes[min(pops[i]):max(pops[i])+1]) for v in sims.variants(as_bytes=True)])
    ps = sstats.PolySIM(sd)
    print(i, ":", ps.tajimasd(),' ',ps.thetaw(),' ',ps.thetapi())

