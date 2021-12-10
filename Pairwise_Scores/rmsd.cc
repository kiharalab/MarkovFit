#include "rmsd.h"

/*
 * This version of the RMSD calculation method assumes that both vectors
 * are the same size and that each position in each array corresponds to
 * the actual alignment of atoms
 */
double calculate_allatom_rmsd(vector<atom>& S, vector<atom>& T)
{
    double rmsd = 0.;
    int N = 0;
    for(size_t atom_index = 0; atom_index < S.size(); atom_index++)
    {
        rmsd += get_squared_distance(S[atom_index].axyz, T[atom_index].axyz);
        N++;
    }
    rmsd /= (double) N;
    rmsd = sqrt(rmsd);
    return rmsd;
}
