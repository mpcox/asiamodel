# example using msprime and libsequence
# https://gist.github.com/molpopgen/519c3f266edaadcc838db278e93fe08d
# code fixed (website version does not run; this does on 24 January 2018)

import msprime as msp
import libsequence.polytable as pt
import libsequence.summstats as sstats
import libsequence.windows as windows

# simple version
g = msp.simulate(sample_size = 10,Ne=1e6, recombination_rate=1e-8,mutation_rate=1e-8,length=1e4)

sd = pt.SimData([(v.position, v.genotypes) for v in g.variants(as_bytes=True)])

ps = sstats.PolySIM(sd)

print(ps.tajimasd(),' ',ps.thetaw(),' ',ps.thetapi())

#Do a "sliding window" analysis, 1 kb at a time, non-overlapping

w = windows.Windows(sd,window_size=1e3,step_len=1e3,starting_pos=0.,ending_pos=1e4)

print(len(w))
for i in w:
    pswi=sstats.PolySIM(i)
    print(pswi.tajimasd(),' ',pswi.thetaw(),' ',pswi.thetapi())



# more complex original version

#Create a list of tuples
#of the format (position, data),
#where data is a 0/1-encoded string
#You can think of this as an "ms block"
#rotated 90 degrees.
#The return value can be input into pylibseq
#to calcualte summary stats
def mksample(g):
    rv = []
    for tree in g.trees():
        muts = list(tree.mutations())
        if len(muts)>0:
            for mi in muts:
                data = ['0']*tree.get_sample_size()
                leaves = tree.leaves(mi.node)
                for li in leaves:
                    data[li]='1'
                rv.append((mi.position,''.join(data)))
    #rv is convientently sorted in order of increasing
    #mutation position, so we don't have to bother with that.
    return rv

g = msp.simulate(sample_size = 10,Ne=1e6, recombination_rate=1e-8,mutation_rate=1e-8,length=1e4)


#sd = pt.SimData([(v.position, v.genotypes) for v in g.variants(as_bytes=True)])


s = mksample(g)

sd = pt.SimData(s)




ps = sstats.PolySIM(sd)

print(ps.tajimasd(),' ',ps.thetaw(),' ',ps.thetapi())

#Do a "sliding window" analysis, 1 kb at a time, non-overlapping

w = windows.Windows(sd,window_size=1e3,step_len=1e3,starting_pos=0.,ending_pos=1e4)

print(len(w))
for i in w:
    pswi=sstats.PolySIM(i)
    print(pswi.tajimasd(),' ',pswi.thetaw(),' ',pswi.thetapi())

