#ifndef _POINT_TRANSFORMATION_H_
#define _POINT_TRANSFORMATION_H_

#include <string>
#include <sstream>
#include <vector>

using std::string;
using std::vector;

// Provides a way of storing a rotation/translation and then applying it to a point, by calling the
// transform method, which returns the x,y,z coordinates after the transformation is applied.
// It also provides a means to undo the transformation

// BIG NOTE: All transformations are in degrees (not radians)
// The class performs the appropriate adjustments internally when libraries require radians
class point_transformation
{
private:
    // Angles that describe the rotation
    // Alpha represents yaw (i.e. z-axis rotation)
    // Beta is the pitch (i.e. y-axis)
    // Gamma is the roll (i.e. x-axis)
    //double alpha, beta, gamma;
    // The "forward" rotation matrix represents a sequence of roll, pitch and yaw rotations, described by the angles provided
    // It's constructuted by multiplying Rz(alpha) * Ry(beta) * Rx(gamma)
    // Where each of the R? matrices are the standard rotation matrices along the corresponding axis
    // Following transformation matrix standards, a new vector
    // v is transformed by multiplying Rv (i.e the matrix multiplies on the left)
    //double rotation_matrix[3][3];
    // Returning the points to their original position is a typical operation in our code.
    // Given that matrix multiplication is not commutative, we precalculate the inverse transformation
    // (i.e. Rx(-gamma) * Ry(-beta) * Rz(-alpha))
    double inverse_rotation_matrix[3][3];
    // After applying rotations, points are translated by these x,y,z magnitudes
    //double translation[3];
    // Auxiliary method to multiply the coordinates by one of the two
    // matrices. gsl could be used but since it's a very specific matrix
    // vector multiplication it's preferred not to include gsl
    void multiply_by_matrix(double matrix[3][3],
                            double in_x, double in_y, double in_z,
                            double* out_x, double* out_y, double* out_z);
public:
    double alpha, beta, gamma;
    double rotation_matrix[3][3];
    double translation[3];
    
    point_transformation() {}
    point_transformation(double alpha_z, double beta_y, double gamma_x,
                         double translate_x, double translate_y,
                         double translate_z);
    // Use this constructor if the rotation matrix is already known, not
    // constructed from the angles given.
    point_transformation(double matrix[3][3], double translate_x,
                         double translate_y, double translate_z);
    point_transformation& operator=(point_transformation other);
    // Accessor properties
    inline double get_alpha() { return alpha; }
    inline double get_beta() { return beta; }
    inline double get_gamma() { return gamma; }
    inline double get_tx() { return translation[0]; }
    inline double get_ty() { return translation[1]; }
    inline double get_tz() { return translation[2]; }
    // Populates both forward and inverse rotation matrices
    void set_rotation(double alpha_z, double beta_y, double gamma_x);
    // Alternate when the matrix is already known.
    void set_rotation(double matrix[3][3]);
    // Simply overwrites the 3D translation vector
    void set_translation(double translate_x, double translate_y, double translate_z);
    // Applies rotation_matrix*[in_x, in_y, in_z] + translation, it first rotates and then translated
    void transform(double in_x, double in_y, double in_z, double* out_x, double* out_y, double* out_z);
    // Inverse process:
    // inverse_rotation_matrix * ([in_x, in_y, in_z] - translation)
    // Note that inverse rotation matrix was constructed using
    // the negative values of the angles provided initially
    void invert_transform(double in_x, double in_y, double in_z, double* out_x, double* out_y, double* out_z);
    // String representation used for printing.
    string to_string();
};

// Used when several transformations are applied to a set of points in sequence.
// e.g. if we have point_transformation t1 and t2, we can add them to
// an instance of this class as follows:
//
// point_transformation_sequence seq;
// seq.add(t1);
// seq.add(t2);
//
// seq.transform(a point) will first apply t1.transform and then t2.transform
// seq.invert_transform(a point) will conversely do t2.invert_trasnform first
// and t1.invert_transform next
class point_transformation_sequence
{
private:
public:
    vector<point_transformation> transformations;
    point_transformation_sequence() {}
    void add(point_transformation transformation);
    void add(point_transformation_sequence sequence);
    // To avoid having the caller create a point_transformation instance if it
    // already has variables for the separate parameters. The 3 methods call
    // add(point_transformation) internally
    void add(double alpha_z, double beta_y, double gamma_x, double translate_x, double translate_y, double translate_z);
    void add_translation(double translate_x, double translate_y, double translate_z);
    void add_rotation(double alpha_z, double beta_y, double gamma_x);
    void clear();
    void transform(double in_x, double in_y, double in_z, double* out_x, double* out_y, double* out_z);
    void invert_transform(double in_x, double in_y, double in_z, double* out_x, double* out_y, double* out_z);
};

// Point transformations can be used to sample the 6D space
// of rotations and translations. Provided the range for each of the
// 3 angles and 3 translational variables, instance of this class
// work as an iterator (in spirit) that replaces 6 nested for loops that
// would be used to generate all combinations of angles and translations.
//
// The typical usage would be to initialize a point_transformation_generator
// with the ranges, followed by calls to next, next, ..., next, where each
// call to next returns a point_transformation instance initialized with the
// current set of angles and translations. A while loop should test for
// is_finished() and continue while this method returns false;
class point_transformation_generator
{
public:
    // All ranges go from xxx_center - xxx_plus_minus to xxx_center
    // to xxx_plus_minus. Each of the variables is updated at a xxx_step rate.
    point_transformation_generator(double alpha_z_rotation_center,
                                   double beta_y_rotation_center, double gamma_x_rotation_center,
                                   double alpha_z_plus_minus, double beta_y_plus_minus,
                                   double gamma_x_plus_minus, double tx_center, double ty_center,
                                   double tz_center, double tx_plus_minus, double ty_plus_minus,
                                   double tz_plus_minus, double rotation_step, double translation_step);
    // When each of the angles and translations are >= center + plus_minus
    bool is_finished();
    // Take the current transformation_parameters, create a point transformation and return it.
    // It will also modify transformation_parameters to have them ready for the next call.
    // Thus this method always computes the transformation that will come next before returning.
    point_transformation next();
    // Returns the total number of transformations that will be returned by next before we finish sampling the space.
    // TODO
    //  unsigned long get_total_transformations();
private:
    // Index from 0 to 5 (0-2 alpha,beta,gamma and 3-5 x,y,z) that indicates which of the six degrees of freedom needs to be updated next
    // As an internal state to indicate when it's finished, the value is set to -1
    int next_index_to_update;
    // The values to be returned as a point_transformation when next is called
    // The order of the values is alpha, beta, gamma, x, y, z
    double transformation_parameters[6];
    // xxx_center - xxx_plus_minus
    double min_values[6];
    // xxx_center + xxx_plus_minus
    double max_values[6];
    // 0-2 -> rotation_step, 3-5 -> translation_step
    double update_steps[6];
};

#endif
