#!/bin/bash

chmod +x *
chmod -R 777 ./*

cd ./output

#Inputs
main_map='../Data/original.mrc'
main_pdb='../Data/1cs4.pdb'
subunit1_map='../Data/A.mrc'
subunit2_map='../Data/B.mrc'
subunit3_map='../Data/C.mrc'

subunit1_pdb='../Data/A.pdb'
subunit2_pdb='../Data/B.pdb'
subunit3_pdb='../Data/C.pdb'

main_contour_level=33
sub_contour_level=33
no_processes=4
voxel_space=1

map_type=2 #2: simulated, 1: experimental
no_clashes=400

#Code
#FFT Search
echo "Search: A"
nohup ../FFT_Search/EMVEC_FIT_PowerFit -a $main_map -b $subunit1_map -t $main_contour_level -T $sub_contour_level -c $no_processes -P true -M 2 -s $voxel_space -p $map_type > A_FFT_search.txt 2>&1;
echo "Search: B"
nohup ../FFT_Search/EMVEC_FIT_PowerFit -a $main_map -b $subunit2_map -t $main_contour_level -T $sub_contour_level -c $no_processes -P true -M 2 -s $voxel_space -p $map_type > B_FFT_search.txt 2>&1;
echo "Search: C"
nohup ../FFT_Search/EMVEC_FIT_PowerFit -a $main_map -b $subunit3_map -t $main_contour_level -T $sub_contour_level -c $no_processes -P true -M 2 -s $voxel_space -p $map_type > C_FFT_search.txt 2>&1;

#Process and Cluster Results
echo "Process: A"
../Handle/handle --dist-threshold 5 --min-dist 8 --correct-x 45.1 --correct-y -14 --correct-z 62.7 --i A_FFT_search.txt --o A_FFT_search_result.txt;
echo "Process: B"
../Handle/handle --dist-threshold 5 --min-dist 8 --correct-x 54.8 --correct-y -16 --correct-z 40.6 --i B_FFT_search.txt --o B_FFT_search_result.txt;
echo "Process: C"
../Handle/handle --dist-threshold 5 --min-dist 8 --correct-x 18.8 --correct-y 4 --correct-z 22.4 --i C_FFT_search.txt --o C_FFT_search_result.txt;

#Compute Pairwise Scores
echo "Pairwise Scores"
nohup mpirun -np $no_processes ../Pairwise_Scores/pairwise_scores_mpi --input-pdb $subunit1_pdb --input-pdb $subunit2_pdb --input-pdb $subunit3_pdb --labels A,B,C --transforms-file A_FFT_search_result.txt --transforms-file B_FFT_search_result.txt --transforms-file C_FFT_search_result.txt  --calpha $main_pdb --transforms-num 100;

#MRF
echo "MRF"
R --slave --vanilla --file=../MRF/mrf.R --args "operation='map'" "singletons=c('A.mrf','B.mrf','C.mrf')" "pairwise=c('A-B.mrf','A-C.mrf','B-C.mrf')" "weights=potential.collection.weights(CC=0.5,Overlap=0.9,PhysicsScore=1,no_clashes=0.8)" "output.prefix='mrf'"

#MaxHeap
echo "MaxHeap"
python3 ../MaxHeap/max-heap.py --mrf-file mrf_top100.txt --clash-threshold $no_clashes -dir ../Data/;

echo "Done!"
$SHELL
