#!/bin/csh -f
pdbset xyzin 5xlr.pdb  xyzout 5xlr_vec.pdb << eof-1
rotate matrix -0.146 0.989 0.000 0.989 0.146 0.000 0.000 0.000 -1.000
shift -160.207 -217.858 174.080
eof-1
