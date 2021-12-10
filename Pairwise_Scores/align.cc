// align.cc
//
// Find rmse of rigid body alignment of target and query, 
// matching residues found by clique analysis
//
// Reference:
//   S. Umeyama, 
//   Least-squares estimation of transformation parameters 
//   between two point patterns,
//   IEEE Transactions on Pattern Analysis and Machine Intelligence,
//   Vol 13, No 4, pp 376-380, 1991.
//
// dw, 17/2/4

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <float.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

using namespace std;


#include "align.h"
#include "pdb.h"


bool debug = 0;


void print_gsl_vector(ostream& os, const gsl_vector *v)
{
    for (size_t i=0; i<v->size-1; i++)
    {
        os << gsl_vector_get(v,i) << " ";
    }
    
    os << gsl_vector_get(v, v->size-1) << endl;
}


void print_gsl_matrix(ostream& os, const gsl_matrix *v)
{
    for(size_t i=0; i<v->size1; i++)
    {
        for(size_t j=0; j<v->size2-1; j++)
        {
            os << gsl_matrix_get(v,i,j) << ", ";
        }
        os << gsl_matrix_get(v,i,v->size2-1) << endl;
    }
}

point_transformation gsl_transformation_to_point_transformation(
                                                                gsl_matrix* R, gsl_vector* t) {
    double rotation[3][3];
    
    for (size_t row = 0; row < 3; row++) {
        for (size_t column = 0; column < 3; column++) {
            rotation[row][column] = gsl_matrix_get(R, row, column);
        }
    }
    
    point_transformation transformation(rotation,
                                        gsl_vector_get(t, 0), gsl_vector_get(t, 1), gsl_vector_get(t, 2));
    return transformation;
}


// find rigid body transformation z=Rx+t minimizing 
// least squares error ||y-z||**2
void align(gsl_matrix_const_view x, gsl_matrix_const_view y, gsl_matrix *R, gsl_vector *t, double *rmse)
{
    
    const size_t n1 = x.matrix.size1;
    const size_t n2 = x.matrix.size2;
    
    gsl_vector *xmu = gsl_vector_alloc(n2);   // means of columns of x
    gsl_vector *ymu = gsl_vector_alloc(n2);   // means of columns of y
    gsl_matrix *X = gsl_matrix_alloc(n1,n2);  // mean-centred x
    gsl_matrix *Y = gsl_matrix_alloc(n1,n2);  // mean-centred y
    
    // for LU decomposition
    int sign;
    gsl_permutation *p = gsl_permutation_alloc(n2);
    gsl_matrix *LU = gsl_matrix_alloc(n2, n2);
    
    // for SVD
    gsl_vector *work = gsl_vector_alloc(n2);
    gsl_matrix *U = gsl_matrix_alloc(n2,n2);
    gsl_matrix *V = gsl_matrix_alloc(n2,n2);
    gsl_vector *d = gsl_vector_alloc(n2);
    
    // mean-centre the data
    for(size_t j=0; j<n2; j++)
    {
        double mu = 0, nu = 0;
        for(size_t i=0; i<n1; i++)
        {
            mu += gsl_matrix_get(&x.matrix,i,j);
            nu += gsl_matrix_get(&y.matrix,i,j);
        }
        gsl_vector_set(xmu,j,mu/n1);
        gsl_vector_set(ymu,j,nu/n1);
    }
    
    if(debug)
    {
        cerr << "xmu: ";
        print_gsl_vector(cerr,xmu);
        cerr << "ymu: ";
        print_gsl_vector(cerr,ymu);
    }
    
    for(size_t i=0; i<n1; i++)
    {
        for(size_t j=0; j<n2; j++)
        {
            gsl_matrix_set(X,i,j,gsl_matrix_get(&x.matrix,i,j)-gsl_vector_get(xmu,j));
            gsl_matrix_set(Y,i,j,gsl_matrix_get(&y.matrix,i,j)-gsl_vector_get(ymu,j));
        }
    }
    
    // form covariance matrix U = Y'X
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0/n1, Y, X, 0.0, U);
    
    // calculate det(covariance matrix) from LU decomposition
    gsl_matrix_memcpy(LU, U);
    gsl_linalg_LU_decomp(LU, p, &sign);
    double det = gsl_linalg_LU_det(LU, sign);
    
    if(debug)
    {
        cerr << "covariance:\n";
        print_gsl_matrix(cerr,U);
        cerr << "det(covariance): " << det << "\n";
    }
    
    // singular value decomposition U -> U*diag(d)*V'
    gsl_linalg_SV_decomp(U, V, d, work) ;
    
    if(debug)
    {
        cerr << "d: ";
        print_gsl_vector(cerr,d);
        cerr << "n2: " << n2 << endl;
    }
    
    // find rank of covariance matrix
    size_t rank;
    
    if(gsl_vector_get(d, n2-1) >= FLT_EPSILON)
    {
        rank = n2;
    }
    else if(gsl_vector_get(d, n2-2) >= FLT_EPSILON)
    {
        rank = n2-1;
    }
    else
    {
        *rmse = -1;
        return;
        //cerr << "align: two zero singular values in svd,\n";
        //exit(EXIT_FAILURE);
    }
    
    
    // calculate detU and detV for rank deficient case
    double detU = 0., detV = 0.;
    
    if(rank == n2-1)
    {
        gsl_matrix_memcpy(LU, U);
        gsl_linalg_LU_decomp(LU, p, &sign);
        detU = gsl_linalg_LU_det(LU, sign);
        
        gsl_matrix_memcpy(LU, V);
        gsl_linalg_LU_decomp(LU, p, &sign);
        detV = gsl_linalg_LU_det(LU, sign);
        
        if(debug)
        {
            cerr << "detU: " << detU << "\n";
            cerr << "detV: " << detV << "\n";
        }
    }
    
    if((rank == n2 && det >= 0) || (rank == n2-1 && detU*detV > 0))
    {
        // R = UV'
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, U, V, 0.0, R);
    }
    else
    {
        // R = USV', S = diag(1,...,1,-1)
        gsl_matrix *S = gsl_matrix_alloc(n2,n2);
        gsl_matrix *T = gsl_matrix_alloc(n2,n2);
        gsl_matrix_set_identity(S);
        gsl_matrix_set(S, n2-1, n2-1, -1.0);
        gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, S, V, 0.0, T);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, T, 0.0, R);
        gsl_matrix_free(S);
        gsl_matrix_free(T);
    }
    
    if(debug)
    {
        cout << "R:\n";
        print_gsl_matrix(cout,R);
    }
    
    // variances
    double sx2 = 0, sy2 = 0;
    
    for(size_t i=0; i<n1; i++)
    {
        for(size_t j=0; j<n2; j++)
            sx2 += gsl_matrix_get(X,i,j) * gsl_matrix_get(X,i,j);
    }
    
    for(size_t i=0; i<n1; i++)
        for(size_t j=0; j<n2; j++)
            sy2 += gsl_matrix_get(Y,i,j) * gsl_matrix_get(Y,i,j);
    
    sx2 /= n1;
    sy2 /= n1;
    
    double trace = 0;
    
    for(size_t i=0; i<n2-1; i++)
        trace += gsl_vector_get(d,i);
    
    if((rank == n2 && det >= 0) || (rank == n2-1 && detU*detV > 0))
    {
        trace += gsl_vector_get(d,n2-1);
    }
    else
    {
        trace -= gsl_vector_get(d,n2-1);
    }
    
    // use next line to include a scaling factor c (z=cRx+t)
    //double c = trace / sx2;
    double c = 1.0;
    
    if(debug)
        cerr << "c: " << c << endl;
    
    // t = -cRxmean + ymean
    gsl_vector_memcpy(t, ymu);
    gsl_blas_dgemv(CblasNoTrans, -c, R, xmu, 1.0, t);
    
    if(debug)
    {
        cout << "t: ";
        print_gsl_vector(cout, t);
    }
    
    // for exact matches rounding errors can produce mse<0
    double mse = sy2 + c*c*sx2 - 2*c*trace;
    *rmse = (mse > 0) ? sqrt(mse) : 0;
    
    // release storage
    gsl_vector_free(xmu);
    gsl_vector_free(ymu);
    gsl_matrix_free(X);
    gsl_matrix_free(Y);
    gsl_matrix_free(LU);
    gsl_vector_free(work);
    gsl_matrix_free(U);
    gsl_matrix_free(V);
    gsl_vector_free(d);
    gsl_permutation_free(p);
}

double align_all_and_transform(vector<atom>* input, vector<atom>* target,
                               vector<atom>* input_to_transform,
                               vector<atom>* transformed_output,
                               point_transformation* optimal_transformation=NULL)
{
    if(input->size() != target->size()) {
        return -1;
    }
    // Create copies needed by GSL
    double* input_coords = new double[input->size() * 3];
    double* target_coords = new double[target->size() * 3];
    
    for(size_t coord_index = 0; coord_index < input->size(); coord_index++) {
        for(size_t xyz = 0; xyz < 3; xyz++) {
            input_coords[coord_index * 3 + xyz] = (*input)[coord_index].axyz[xyz];
            target_coords[coord_index * 3 + xyz] = (*target)[coord_index].axyz[xyz];
        }
    }
    // The copies are fed into a view
    gsl_matrix_const_view input_view = gsl_matrix_const_view_array(input_coords, input->size(), 3);
    gsl_matrix_const_view target_view = gsl_matrix_const_view_array(target_coords, target->size(), 3);
    
    // Calculate the geometric transformation and get the RMSD
    gsl_matrix *R = gsl_matrix_alloc(3,3);
    gsl_vector *t = gsl_vector_alloc(3);
    double rmse;
    align(input_view, target_view, R, t, &rmse);
    
    if(rmse != -1) {
        if(debug) {
            cerr  << "Transformation matrix:\n"
            << gsl_matrix_get(R,0,0) << "\t" << gsl_matrix_get(R,1,0)
            << "\t" << gsl_matrix_get(R,2,0) << "\n"
            << gsl_matrix_get(R,0,1) << "\t" << gsl_matrix_get(R,1,1)
            << "\t" << gsl_matrix_get(R,2,1) << "\n"
            << gsl_matrix_get(R,0,2) << "\t" << gsl_matrix_get(R,1,2)
            << "\t" << gsl_matrix_get(R,2,2) << "\n"
            << "\nTranslation:\n" << gsl_vector_get(t,0) << "\t"
            << gsl_vector_get(t,1) << "\t" << gsl_vector_get(t,2)
            << endl << endl;
            cerr << "Starting transformation of " << transformed_output->size()
            << " atoms\n";
        }
        
        // If the alignment was successful then create the transformed coordinates
        for(size_t transformed_index = 0; transformed_index < input_to_transform->size(); transformed_index++) {
            atom transformed_atom = (*input_to_transform)[transformed_index];
            if(debug) {
                cerr << "Transforming " << transformed_atom.residue << " " << transformed_atom.rnum << endl;
            }
            
            double new_coords[3]; //x,y,z
            for(size_t axis = 0; axis < 3; axis++)
            {
                new_coords[axis] =
                transformed_atom.axyz[0] * gsl_matrix_get(R,axis,0) +
                transformed_atom.axyz[1] * gsl_matrix_get(R,axis,1) +
                transformed_atom.axyz[2] * gsl_matrix_get(R,axis,2) +
                gsl_vector_get(t,axis);
            }
            // Done twice because we use all 3 coordinates in each axis computation
            for(size_t axis = 0; axis < 3; axis++) {
                transformed_atom.axyz[axis] = new_coords[axis];
            }
            transformed_output->push_back(transformed_atom);
        }
    }
    
    // The caller can provide a pointer where a copy of the transformation
    // is stored.
    if (optimal_transformation) {
        *optimal_transformation = gsl_transformation_to_point_transformation(R, t);
    }
    
    gsl_matrix_free(R);
    gsl_vector_free(t);
    delete [] input_coords;
    delete [] target_coords;
    
    return rmse;
}
