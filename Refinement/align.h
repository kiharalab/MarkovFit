// align.h
//
// prototypes for align.cc
//
// dw, 12/12/3

#include <gsl/gsl_matrix.h>

#include "pdb.h"
#include "point_transformation.h"
#include <map>

// find rigid body transformation z=Rx+t minimizing least squares error ||y-z||**2
void align(gsl_matrix_const_view , gsl_matrix_const_view, gsl_matrix *, gsl_vector *, double *);
// Not only return the RMSD value of the alignment but it also returns a new
// vector with the modified x,y,z coordinates that make input align on target.
// As the method name states ALL atoms are aligned. It's expected that both
// collections are the same size.
// optimal_transformation is an optional output parameter that will copy
// the rotation matrix and translation used to obtain the alignment
// Note: The reason why there are 2 inputs is because the first one can be
// just c-alphas, for instance, and the latter one can be an all atom
// representation that needs to be transformed.
double align_all_and_transform(vector<atom>* input, vector<atom>* target,
    vector<atom>* input_to_transform, vector<atom>* transformed_output,
    point_transformation* optimal_transformation);
