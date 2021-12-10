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
##### Input:
##### Output:

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

