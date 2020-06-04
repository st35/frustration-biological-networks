import warnings
warnings.filterwarnings("ignore")

Sample_Names = []
Sample_Type = []
with open('dump/Sample_Description.log', 'r') as f:
    count = 0
    for line in f:
        l = line.strip().split('\t')
        l = l[1:]
        if count == 0:
            for val in l:
                Sample_Names.append(val.strip())
            count += 1
            continue
        if count == 1:
            for val in l:
                Sample_Type.append(val.strip())
            count += 1
Normal_Samples = []
for i in range(len(Sample_Names)):
    if '"tissue type: Non-tumor"' in Sample_Type[i]:
        Normal_Samples.append(Sample_Names[i][1:-1])

import pandas as pd

D = pd.read_csv('dump/Terenuma_Data.log', sep = '\t')
Normal_Data = D[Normal_Samples]
avail_genes = list(D.index)

MYC_genes = []
with open('dump/MYC_genes.log', 'r') as f:
    for line in f:
        MYC_genes.append(line.strip())

avail_MYC_genes = [gene for gene in MYC_genes if gene in avail_genes]

genes_to_use = set()
TFs = set()
with open('dump/TRRUST_Human_Data.log', 'r') as f:
    for line in f:
        l = line.strip().split('\t')
        source = l[0].strip()
        target = l[1].strip()
        TFs.add(source)
        if target in avail_MYC_genes:
            if source in avail_genes:
                genes_to_use.add(source)
                genes_to_use.add(target)

import numpy as np
from arboreto.algo import grnboost2, genie3
from sklearn.decomposition import PCA

pca = PCA(n_components = 1)

if __name__ == '__main__':
    pcafile = open('pc_fraction_explained.log', 'w')
    data_to_use = Normal_Data.loc[list(genes_to_use), :].T
    pca.fit(data_to_use)
    pca_explained = pca.explained_variance_ratio_[0]
    pcafile.write(str(pca_explained) + '\n')
    TFs_to_use = [gene for gene in genes_to_use if gene in TFs]
    network = grnboost2(expression_data = data_to_use, tf_names = list(TFs_to_use))
    print(network.head())
    network.to_csv('networkfiles/biologicalnetwork.log', sep = '\t', index = False, header = False)

    original_data = Normal_Data.loc[list(genes_to_use), :]
    for i in range(100):
        print(i)
        data_to_use = original_data.copy()
        for j in range(len(Normal_Samples)):
            l = list(original_data.iloc[:, j])
            l = list(np.random.permutation(l))
            data_to_use[Normal_Samples[j]] = l
        print(data_to_use.shape)
        data_to_use = data_to_use.T
        pca.fit(data_to_use)
        pca_explained = pca.explained_variance_ratio_[0]
        pcafile.write(str(pca_explained) + '\n')
        TFs_to_use = [gene for gene in genes_to_use if gene in TFs]
        network = grnboost2(expression_data = data_to_use, tf_names = list(TFs_to_use))
        print(network.head())
        network.to_csv('networkfiles/networkrandomdata_' + str(i + 1) + '.log', sep = '\t', index = False, header = False)
    pcafile.close()
