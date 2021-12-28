# MarkovFit
MarkovFit is A Protein Structural Modeling Method for Electron Microscopy Maps Using Markov Random Field.

Copyright (C) 2022 Eman Alnabati, Genki Terashi, Juan Esquivel-Rodriguez, Daisuke Kihara, and Purdue University.

License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)

Contact: Daisuke Kihara (dkihara@purdue.edu)

Cite:
## Outline
- Pre-required software
    - Compile Source Code
- Compile Source Code
- Project Steps
    - Generate simulated maps of subunits using EMAN2 package
    - Run FFT Search
    - Handle and Cluster FFT Search Results
    - Compute Pairwise Scores of Pairs of Subunits
    - Generate MRF graph and Apply Belief Propagation
    - Extract top10 Final Structures using MaxHeap Tree
    - Physics-based Refinement
- Run an example
    

## Pre-required software
- Python 3 : https://www.python.org/downloads/
- EMAN2 : https://blake.bcm.edu/emanwiki/EMAN2/Install/
- Chimera : https://www.cgl.ucsf.edu/chimera/download.html
- FFTW: http://www.fftw.org/download.html
- GCC Compiler

### Compile Source Code:
- FFT_Search Folder:
make EMVEC_FIT_PowerFit

- Handle Folder:
g++ handle_result.cc -o handle

- Pairwise_Scores Folder:
make pairwise_scores_mpi

- Refinement Folder:
make refine_mpi

## Project Steps:
### Generate simulated maps of subunits using EMAN2 package:
    e2pdb2mrc.py --apix=voxel_spacing --res=resolution input_pdb_file output_mrc_file

For experimental map fitting: use voxel spacing and resoultion of that map. 

### Run FFT Search:
##### Command:
```sh
./FFT_Search/EMVEC_FIT_PowerFit -a main_map -b subunit1_mrc_map -t main_map_contour_level -T subunit_map_contour_level -c no_processes -P true -M 2 -s voxel_space -p map_type > output_file
```

##### Input:
    -a:          Main map
    -b:          Subunit map
    -t [float]:  Threshold of density main_map def=0.000
    -T [float]:  Threshold of density sub_map def=0.000
    -c [int  ]:  Number of cores for threads def=2
    -g [float]:  Bandwidth of the gaussian filter
                 def=16.0, sigma = 0.5*[float]
    -s [float]:  Sampling grid space def=7.0
    -M [int]:    Sampling Angle interval Mode 1-3 def=2
                 1: 20.83 degree,   648 samples
                 2: 10.07 degree, 7,416 samples
                 3: 4.71 degree, 70,728 samples
    -C:          Cross Correlation Coefficient and Overlap Mode 
                 Using normalized density value by Gaussian Filter
    -P:          Pearson Correlation Coefficient and Overlap Mode 
                 Using normalized density value by Gaussian Filter and average density
    -p:          Map type: 1 for experimental, 2 for simulated def=1 
           
##### Output:
    output_file contains the different transformations applied to subunit map along with goodness-of-fit scores. 

### Handle and Cluster FFT Search Results:
##### Command:
```sh
./Handle/handle --dist-threshold max_dist_to_show --min-dist cluster_dist_thrshold --correct-x center_x --correct-y center_y --correct-z center_z --i input_file --o output_file
```
##### Input:
    --min-dist:         Distance used for clustering search results 8 def=8
    --i:                Input file containing FFT search results
    --o:                Output file to store results after sorting and clustering
    optional:
    To print search results and their distance to reference structur after sorting and clustering
    --correct-x:        X value of center of reference/native structure
    --correct-y:        Y value of center of reference/native structure 
    --correct-z:        Z value of center of reference/native structure 
    --dist-threshold:   Show only results with distance < threshold to native structure

##### Output:
    output_file containing sorted and clustered results 
    output_file_b4_clstering containing sorted results before clustering

### Compute Pairwise Scores of Pairs of Subunits
##### Command:
```sh
 mpirun -np no_processes ./Pairwise_Scores/pairwise_scores_mpi --input-pdb subunit1_pdb ... --input-pdb subunitN_pdb --labels A,B,C --transforms-file sub1_processed_search_result_file --transforms-file subN_processed_search_result_file  --calpha main_pdb --transforms-num no_results_per_subunit
```
##### Input:
    -np:                    No of processes
    --input-pdb:            Subunit pdb files
    --labels:               List of subunits IDs separated by comma
    --transforms-file:      Subunit processed search result file 
    --calpha main_pdb:      Native structure pdb file
    --output-prefix:        Prefix for output files
    --transforms-num:       No results to be considered for each subunit

##### Output:
    One .mrf file for each subunit containing RMSD and goodness-of-fit scores for each transformation
    One .mrf file for each pair of subunits containing RMSD and pairwise scores for each of their transformations
    
### Generate MRF graph and Apply Belief Propagation
##### Command:
```sh
R --slave --vanilla --file=./MRF/mrf.R --args "operation='map'" "singletons=c('sub1.mrf',...,'subN.mrf')" "pairwise=c('sub1-sub2.mrf','sub1-sub3.mrf',...,'sub(N-1)-subN.mrf')" "weights=potential.collection.weights(CC=0.5,Overlap=0.9,PhysicsScore=1,no_clashes=0.8)" "output.prefix='mrf'"
```
##### Input:
    singletons:         List of .mrf files of all subunits
    pairwise:           List of .mrf files of all pairs of subunits
    weights:            weights to be used for the goodness-of-fit scores and pairwise scores  
    output.prefix:      Prefix of the output file

##### Output:
    prefix_top100.txt file containing the transformations of subunits and pairs of subunits sorted bt their final beliefs computed by belief propagation algorithm.

### Extract top10 Final Structures using MaxHeap Tree:
##### Command:
```sh
python3 ./MaxHeap/max-heap.py --mrf-file prefix_top100.txt --clash-threshold no_clashes -dir pdb_dir
```
##### Input:
    --mrf-file:             File containing the results of MRF and belief propagation 
    --clash-threshold:      Maximum number of acceptable clashes between pairs of subunits
    -dir:                   Directory of the native structure pdb files to generate the transformed structures   

##### Output:
    PDB files of the top 10 structures of each subuint (subID_decoy_MaxHeap_structure#.pdb)
    File for each subunit containing subunit transformations of the top 10 structures (subID_clashesThreshold_MaxHeap.txt)
    
### Physics-based Refinement
##### Command:
```sh
 mpirun -np no_processes ./Refinement/refine_mpi --input-pdb subunit1_pdb ... --input-pdb subunitN_pdb --labels A,B,C --transforms-file sub1_ID_400_MaxHeap.txt ... --transforms-file subN_ID_400_MaxHeap.txt --calpha main_pdb
```
##### Input:
    -np:                    No of processes max=10
    --input-pdb:            Subunit pdb files 
    --labels:               List of subunits IDs separated by comma
    --transforms-file:      Subunit MaxHeap result file
    --calpha:               Native structure pdb file
    --output-prefix:        Prefix for output files


##### Output:
        PDB files of the refined top 10 structures of each subuint (subID_phy_refine_decoy_structure#.pdb)

## Run an example:
##### Command:
```sh
 ./ run.sh
```
run.sh file contains all the commands to geenrate top 10 refined structures of the target in the output folder.
Expected_Output folder contains all output files of the example target. 
- subunitID_decoy_MaxHeap_#.pdb: PDB structures before refinement.
- subunitID_phy_refine_decoy_#.pdb : PDB structures after refinement.
