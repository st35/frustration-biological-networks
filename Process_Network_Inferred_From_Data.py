import sys

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
import numpy as np

D = pd.read_csv('dump/Terenuma_Data.log', sep = '\t')
Normal_Data = D[Normal_Samples]

for i in range(1, 101):
    Chosen_Pairs = []
    with open('networkfiles/networkrandomdata_' + str(i) + '.log', 'r') as f:
        count = 0
        for line in f:
            l = line.strip().split('\t')
            Chosen_Pairs.append([])
            Chosen_Pairs[count].append(l[0].strip())
            Chosen_Pairs[count].append(l[1].strip())
            count += 1
            if count == 500:
                break
    Sources = set()
    Targets = set()
    Genes = set()
    Sign = []
    for j in range(len(Chosen_Pairs)):
        Sources.add(Chosen_Pairs[j][0])
        Targets.add(Chosen_Pairs[j][1])
        Genes.add(Chosen_Pairs[j][0])
        Genes.add(Chosen_Pairs[j][1])
        gene1 = list(Normal_Data.loc[Chosen_Pairs[j][0]])
        gene2 = list(Normal_Data.loc[Chosen_Pairs[j][1]])
        c = np.corrcoef(gene1, gene2)
        if c[0, 1] > 0:
            Sign.append(1)
        elif c[0, 1] < 0:
            Sign.append(2)
    To_remove = []
    for gene in Genes:
        if gene not in Sources or gene not in Targets:
            To_remove.append(gene)
    print(str(len(Chosen_Pairs)) + '\t' + str(len(Sign)))
    lineswritten = 0
    with open('inputfiles/networkrandomdata' + '_' + str(i) + '.topo', 'w') as f:
        f.write('Source\tTarget\tType\n')
        for j in range(len(Chosen_Pairs)):
            source = Chosen_Pairs[j][0]
            target = Chosen_Pairs[j][1]
            if source in To_remove or target in To_remove:
                continue
            f.write(source + '\t' + target + '\t' + str(Sign[j]) + '\n')
            lineswritten += 1
    IDS = []
    with open('inputfiles/networkrandomdata' + '_' + str(i) + '.topo', 'r') as f:
        count = 0
        for line in f:
            if count == 0:
                count += 1
                continue
            l = line.strip().split('\t')
            source = l[0].strip()
            target = l[1].strip()
            if source not in IDS:
                IDS.append(source)
            if target not in IDS:
                IDS.append(target)
    with open('inputfiles/networkrandomdata' + '_' + str(i) + '.ids', 'w') as f:
        for j in range(len(IDS)):
            f.write(IDS[j] + '\t' + str(j) + '\n')
    print(len(IDS))
    with open('inputfiles/networkrandomdata' + '_' + str(i) + '.phs', 'w') as f:
        f.write('0\t1\n')
