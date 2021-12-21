import numpy as np
import matplotlib.pyplot as plt

# %matplotlib inline

import bioframe
import pandas as pd
import cooler

import cooltools
import cooltools.eigdecomp
import cooltools.expected
import cooltools.saddle

# from dask.distributed import Client, LocalCluster
# cluster = LocalCluster(n_workers=2, processes=2)
# client = Client(cluster)

import re
import sys

coolpath = sys.argv[1]
sampleName = sys.argv[2]
# coolpath = './100kb_cool_file/MergeHiCPro_allValidPairs.consistent.100000.cool'
c = cooler.Cooler(coolpath)

## define region
regions = [(chrom, 0, c.chromsizes[chrom]) for chrom in c.chromnames]
bins = c.bins()[:]
genecov = bins

# import pandas as pd
## the gene count in 100kb bins
bins_cov = pd.read_csv('/Users/wufan/Documents/work/20190910_HiC_hm_terminal_cellCycle_LD_MergeAll/cool_file/100kb_cool_file/cov.bed',sep='\t',header=None)
genecov['gene_count'] = bins_cov.iloc[:,-4]
genecov['gene_coverage'] = bins_cov.iloc[:,-1]

# Perform eigenvector decomposition in cis, sorting and flipping eigenvectors
# according to their correlation with the number of genes in each bin.
cis_eigs = cooltools.eigdecomp.cooler_cis_eig(
    c,
    genecov,
    regions=None,
    n_eigs=5, ## top 5 eigenvectors
    phasing_track_col='gene_count')


Q_LO = 0.025 # ignore 2.5% of genomic bins with the lowest E1 values
Q_HI = 0.975 # ignore 2.5% of genomic bins with the highest E1 values
N_GROUPS = 38 # divide remaining 95% of the genome into 38 equisized groups, 2.5% each
q_edges = np.linspace(Q_LO, Q_HI, N_GROUPS+1)


# Filter track used for grouping genomic bins based on bins filtered out in Hi-C balancing weights
# Doesn't do anything with eigenvectors from the same Hi-C data (hence commented out here),
# but important for external data, such as ChIP-seq tracks
eig = cooltools.saddle.mask_bad_bins((cis_eigs[1], 'E1'), (c.bins()[:], 'weight'))
group_E1_bounds = cooltools.saddle.quantile(eig['E1'], q_edges)
print(group_E1_bounds)

# Assign the group to each genomic bin according to its E1, i.e. "digitize" E1.
digitized, hist = cooltools.saddle.digitize_track(
    group_E1_bounds,
    track=(eig, 'E1'),
)

# Calculate the decay of contact frequency with distance (i.e. "expected")
# for each chromosome.
expected = cooltools.expected.cis_expected(c, regions, use_dask=True)

# Make a function that returns observed/expected dense matrix of an arbitrary
# region of the Hi-C map.
getmatrix = cooltools.saddle.make_cis_obsexp_fetcher(c, (expected, 'balanced.avg'))


# Compute the saddle plot, i.e. the average observed/expected between genomic
# ins as a function of their digitized E1.
S, C = cooltools.saddle.make_saddle(
    getmatrix,
    group_E1_bounds,
    (digitized, 'E1' + '.d'),
    contact_type='cis')

# interaction_sum : 2D array
#     The matrix of summed interaction probability between two genomic bins
#     given their values of the provided genomic track.
# interaction_count : 2D array
#     The matrix of the number of genomic bin pairs that contributed to the
#     corresponding pixel of ``interaction_sum``.


from matplotlib.backends.backend_pdf import PdfPages

pdf = PdfPages(sampleName + '.pdf')
plt.imshow(
    np.log2(S / C)[1:-1, 1:-1],
    cmap='coolwarm',
    vmin=-1,
    vmax=1,
)
plt.colorbar(label='log2 obs/exp')
print('savefig...')
pdf.savefig()
plt.close()
pdf.close()

# outDF.iloc[0:8,0:8] ## BB
# outDF.iloc[0:8,-8:] ## BA
# outDF.iloc[-8:,0:8] ## AB
# outDF.iloc[-8:,-8:] ## AA

outDF = pd.DataFrame((S / C)[1:-1,1:-1])
outDF.to_csv(sampleName+'log2_OBSEXP.csv',header=False,index=False)

print('Some data:')
print('\n')
print(sampleName+' compartment strength:')
print((np.median(outDF.iloc[0:8,0:8]) + np.median(outDF.iloc[-8:,-8:]))/(np.median(outDF.iloc[0:8,-8:])+np.median(outDF.iloc[-8:,0:8])))

print(sampleName+' AA interaction:')
print(np.log2(np.median(outDF.iloc[-8:,-8:])))
print(sampleName+' BB interaction:')
print(np.log2(np.median(outDF.iloc[0:8,0:8])))
print(sampleName+' AB/BA interaction:')
print(np.log2(np.median(outDF.iloc[-8:,0:8])))





