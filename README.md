# MarkovFit
MarkovFit is A Protein Structural Modeling Method for Electron Microscopy Maps Using Markov Random Field.

Copyright (C) 2022 Eman Alnabati, Genki Terashi, Juan Esquivel-Rodriguez, Daisuke Kihara, and Purdue University.

License: GPL v3 for academic use. (For commercial use, please contact us for different licensing.)

Contact: Daisuke Kihara (dkihara@purdue.edu)

Cite:

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
e2pdb2mrc.py --apix=voxel_spacing --res=resolution input_pdb_file output_mrc_file;

For experimental map fitting: use voxel spacing and resoultion of that map. 

### Run FFT Search:
##### Command:
./FFT_Search/EMVEC_FIT_PowerFit -a main_map -b subunit1_mrc_map -t main_map_contour_level -T subunit_map_contour_level -c no_processes -P true -M 2 -s voxel_space -p map_type > output_file;
##### Input:
-a         : Main map
-b         : Subunit map
-t [float] : Threshold of density main_map def=0.000
-T [float] : Threshold of density sub_map def=0.000
-c [int  ] : Number of cores for threads def=2
-g [float] : Bandwidth of the gaussian filter
             def=16.0, sigma = 0.5*[float]
-s [float] : Sampling grid space def=7.0
-M [int]   : Sampling Angle interval Mode 1-3 def=2
             1: 20.83 degree,   648 samples
             2: 10.07 degree, 7,416 samples
             3: 4.71 degree, 70,728 samples
-C         : Cross Correlation Coefficient and Overlap Mode 
             Using normalized density value by Gaussian Filter
-P         : Pearson Correlation Coefficient and Overlap Mode 
             Using normalized density value by Gaussian Filter and average density
-p         : Map type: 1 for experimental, 2 for simulated def=1 
           
##### Output:
output_file contains the different transformations applied to subunit map along with goodness-of-fit scores. 

### Handle and Cluster FFT Search Results:
##### Command:
##### Input:
##### Output:

### Compute Pairwise Scores of Pairs of Subunits
##### Command:
##### Input:
##### Output:

### Generate MRF graph and Apply Belief Propagation
##### Command:
##### Input:
##### Output:

### Extract top10 Final Structures using MaxHeap Tree:
##### Command:
##### Input:
##### Output:

### Physics-based Refinement
##### Command:
##### Input:
##### Output:

