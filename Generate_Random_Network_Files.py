import sys

R = []
with open('../networkfiles/' + sys.argv[1] + '.topo', 'r') as f:
	count = 0
	for line in f:
		if count == 0:
			count += 1
			continue
		l = line.strip().split('\t')
		R.append([l[0].strip(), l[1].strip(), l[2].strip()])

from math import floor
from math import log
import random
from datetime import datetime
random.seed(datetime.now())

for i in range(750):
	print i
	R_new = []
	for j in range(len(R)):
		R_new.append([R[j][0], R[j][1], R[j][2]])
	N = int(floor(log(1e7)*len(R_new) / 2))
	for j in range(N):
		edge1 = int(floor(random.randint(0, len(R_new) - 1)))
		edge2 = edge1
		while edge2 == edge1:
			edge2 = int(floor(random.randint(0, len(R_new) - 1)))
		target1 = R_new[edge1][1]
		target2 = R_new[edge2][1]
		R_new[edge1][1] = target2
		R_new[edge2][1] = target1
		
	if i == 0:
		with open('inputfiles/' + sys.argv[1] + '_' + str(i) + '.topo', 'w') as f:
			f.write('Source\tTarget\tType\n')
			for j in range(len(R)):
				f.write(R[j][0] + '\t' + R[j][1] + '\t' + R[j][2] + '\n')
	else:
		with open('inputfiles/' + sys.argv[1] + '_' + str(i) + '.topo', 'w') as f:
			f.write('Source\tTarget\tType\n')
			for j in range(len(R_new)):
				f.write(R_new[j][0] + '\t' + R_new[j][1] + '\t' + R_new[j][2] + '\n')
	with open('../networkfiles/' + sys.argv[1] + '.ids', 'r') as f, open('inputfiles/' + sys.argv[1] + '_' + str(i) + '.ids', 'w') as g:
		for line in f:
			g.write(line.strip() + '\n')
	with open('../networkfiles/' + sys.argv[1] + '.phs', 'r') as f, open('inputfiles/' + sys.argv[1] + '_' + str(i) + '.phs', 'w') as g:
		for line in f:
			g.write(line.strip() + '\n')
