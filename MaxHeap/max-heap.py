#!/usr/bin/env python
import sys
import numpy as np
import heapq 
import os
import argparse



# parse the arguments from command line
parser = argparse.ArgumentParser(description = 'Output top 10 mrf structures using MaxHeap tree. python3 max-heap.py --mrf-file mrf_result-file --clash-threshold 400 -dir pdb-files-dir.');

parser.add_argument('-mrf', '--mrf-file', required = True, action = 'store', dest = 'mrf_result', help = 'Required. Name of the mrf result file.')
parser.add_argument('-clashes', '--clash-threshold', required = True, action = 'store',type=int, dest = 'clash_threshold', help = 'Required. Max number of clashes accepted between structure subunits.')
parser.add_argument('-dir', required = True, action = 'store', dest = 'pdb_dir', help = 'Directory where the PDB files are located.')

args = parser.parse_args()
#argparse.FileType('r')

#parser.add_argument('file', type=argparse.FileType('r'))
#args = parser.parse_args()

#print(args.file.readlines())



f = args.mrf_result
with open(args.mrf_result, 'r') as file:
	flines=file.readlines()
clash_threshold  = args.clash_threshold
# print(clash_threshold)
pdb_dir = args.pdb_dir
int(clash_threshold)

lines=[]
data = []
pairwise_data = []

# flines =x.readlines()
# print(flines)
#flines = args.mrf_result.readlines()
#print(flines)
singleton_flag = True

#Read data
for i in flines:
    l = i.split(None,8)
    #print(l[0])
    if l[0][0].isalnum()==True  and "-" in l[0]:
        singleton_flag = False
    if l[0].isalnum()==False:
        if singleton_flag == True:
            l[8] = l[8].rstrip()
            res = [float(idx) for idx in l[8].split()]

            data.append(l[:8]+ [str(res[0])])

        if singleton_flag == False and "-" not in l[0]:
            ll = i.split(None,16)
            ll[16] = ll[16].rstrip()
            pairwise_data.append(ll)

a = np.array(data)
#Find singleton names and counts
(unique, counts) = np.unique(a[:,7], return_counts=True)
no_chains = np.count_nonzero(unique)

#Find pair names and counts
b = np.array(pairwise_data)
b_left = np.array(b[:,1:7],dtype='float')
b_right = np.array(b[:,8:14],dtype='float')

pairs=np.apply_along_axis(lambda d: d[7] + '-' + d[14], 1, b)
(pair_unique, pair_counts) = np.unique(pairs, return_counts=True)
no_pairs = np.count_nonzero(pair_unique)

#Get starting index of every chain
no_beliefs = np.max(counts)
result = no_beliefs == np.min(counts)
beliefs = np.zeros((no_beliefs,no_chains+1))
start_index = np.unique(a[:,7],return_index = True)[1]

#Split Beliefs based on chain
temp = np.split(a[:,8], start_index)
del temp[0]

#Populate beliefs
count = 0
for xi in temp:
    xi=np.float_(xi)
    if len(xi)<no_beliefs:
        xi = np.append( xi,np.zeros(no_beliefs-len(xi)))
    beliefs[:,count] = xi
    count+=1

#Add sum of beliefs coluomn
beliefs[:,no_chains] = np.sum(beliefs,axis = 1)

no_results = 10

#Indicies of of each sum of beleifs
indices = np.zeros((no_beliefs,no_chains+1))
indices[:,no_chains] = beliefs[:,no_chains]

RR = np.empty((0,no_chains))

tenth = np.array(range(0, no_beliefs))

for i in range(no_chains):
    indices[:,i] = tenth

belief = np.multiply(beliefs[:,-1],-1).tolist()

#Create max-heap
heapq.heapify(belief) 

items = [0,1,2,3,4,5,6]

print('Finding conformations:')
#max-heap part
i = 1
while i <= no_results :
    if len(belief) ==0:
        print('\n\nNo Belief!\n\n')
        break
    #print('Result no. ', i)
    #extract max
    max = heapq.heappop(belief)*-1
    #print('max = ',max)
    #find max-belief indices
    result = np.where(indices[:,-1] == max)
    new_indices = np.full(no_chains,-1)
    #check max-belief is found: yes
    if len(result) > 0 and len(result[0]) > 0:
        indices_row = result[0][0]
        #extract and print chain transformations
        chain_temp = np.empty((no_chains, 9),dtype='U25')

        for j in range(no_chains):
            chain_temp[j] = data[int(indices[indices_row,j])+start_index[j]]
            new_indices[j] = indices[indices_row,j]

        #Check for clashes
        pair_temp = np.empty((no_pairs, 2))
        clash_flag = 0
        for j in range(no_pairs):
            left = pair_unique[j][0]
            right = pair_unique[j][2]

            m1 = chain_temp[np.where(chain_temp[:,7]==left),1:7][0][0]
            
            m1 = [float(i) for i in m1]

            m2 = chain_temp[np.where(chain_temp[:,7]==right),1:7][0][0]
            m2 = [float(i) for i in m2]

            pair_temp[j][0] = np.intersect1d(np.intersect1d(np.where(b[:,7] == left)[0],np.where(b[:,14] == right)[0]),
                                         np.intersect1d(np.nonzero(np.all(np.where(b_left == m1,1,0), axis=1))[0],
                                                        np.nonzero(np.all(np.where(b_right == m2,1,0), axis=1))[0]))[0]
                                                        
            pair_temp[j][1] = b[int(pair_temp[j][0]),15]

            # Check for clashes
            if pair_temp[j][1] > clash_threshold:
                clash_flag = 1
                #print(i, ' CLASH ',' indices ',new_indices,'  # clashes =  ', pair_temp[j][1])
                break
                    
        #Conformation doesn't have clashes, print it
        if clash_flag == 0:
            found = np.nonzero(np.all(np.where(RR == new_indices,1,0), axis=1))[0]

            if found.size > 0:
                print(new_indices)
                continue
        
            RR= np.vstack([RR, new_indices])

            print('Result no. ', i)
            print(' Sum of Beliefs = ',max)
            print(' Subunit indices ',new_indices)

            for j in range(no_chains):
                print('chain ', j, ' index ', indices[indices_row,j], ' ', chain_temp[j])
                chain_file = unique[j]+'_'+str(clash_threshold)+'_MaxHeap.txt'

                if i==1:
                    with open(chain_file, 'w') as f:
                        for item in items:
                            f.write("%s\t" % chain_temp[j,item])
                        f.write("\n")
                else:
                    with open(chain_file, 'a') as f:
                        for item in items:
                            f.write("%s\t" % chain_temp[j,item])
                        f.write("\n")
            i= i+1
            print("\n")

        #add new beliefs
        for j in range(no_chains):
            if new_indices[j]+1 < no_beliefs:
                new_indices[j] += 1
                sum_b = 0
                flag = 0
                for w in range(len(indices)):
                    comparison = indices[w,:-1] == new_indices
                    equal_arrays = comparison.all()
                    if equal_arrays:
                        flag = 1

                if flag == 0:
                    for k in range(no_chains):
                        sum_b += beliefs[int(new_indices[k]),k]
                    sum_b *=-1

                    #add new beleif
                    heapq.heappush(belief, sum_b)
                    new_indices = np.append(new_indices, sum_b*-1)
                    indices = np.vstack([indices, new_indices])
                    new_indices = np.delete(new_indices,-1)
                
                new_indices[j] -= 1
f.close()

chains = ','.join(unique)
# unique.with(,)
# print(os. getcwd() )
# print(sys.path[0])
command = 'python3 '+(sys.path[0])+'/transform_pdbs.py -i2 ' + chains + ' -t _400_MaxHeap.txt -output decoy_MaxHeap -odir '+pdb_dir
os.system(command)


with open('result.txt', 'w') as f:
    for item in RR:
        f.write("%s\n" % item)
