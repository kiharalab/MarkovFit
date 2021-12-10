#ifndef _RMSD_H_
#define _RMSD_H_

#include <algorithm>
#include <map>
#include <cmath>
#include <iostream>
#include <vector>
#include "pdb.h"
#include "atom.h"

// Indicates that, when calculating the ligand rmsd between a receptor/ligand pair
// only atoms (in the ligand) closer than 10 angstroms to any receptor atom should
// be considered
#define LIGAND_RMSD_THRESHOLD 10

/*
 * This version of the RMSD calculation method assumes that both vectors
 * are the same size and that each position in each array corresponds to
 * the actual alignment of atoms
 */
double calculate_allatom_rmsd(vector<atom>& S, vector<atom>& T);

#endif
