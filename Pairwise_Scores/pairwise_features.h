#ifndef _PAIRWISE_FEATURES_H_
#define _PAIRWISE_FEATURES_H_

#include <fstream>
#include "soroban_score.h"

using std::ofstream;
using std::string;

class mrf_one_body_features {
public:
    mrf_one_body_features(double alpha_z, double beta_y, double gamma_x, double translate_x, double translate_y,
                          double translate_z, double rmsd, double correlation, double overlap, string label);
    void write(ofstream& data_stream);
    // Tab separated list of column names
    static void write_header(ofstream& data_stream);
private:
    // Identifies the transformations that originated these feature values
    double alpha_z, beta_y, gamma_x, translate_x, translate_y, translate_z;
    // Derived from the difference between the transformed unit and the c-alpha trace.
    // This can later be used to label positive and negative cases.
    double rmsd;
    // Derived from the unit EM map and the complete map.
    double correlation, overlap;
    string label;
};

class mrf_two_body_features {
public:
    mrf_two_body_features(double alpha_z_left, double beta_y_left, double gamma_x_left, double translate_x_left,
                          double translate_y_left, double translate_z_left, double alpha_z_right, double beta_y_right,
                          double gamma_x_right, double translate_x_right, double translate_y_right, double translate_z_right,
                          double rmsd, soroban_score physics_score, int no_clashes, string label_left, string label_right);
    void write(ofstream& data_stream);
    // Tab separated list of column names
    static void write_header(ofstream& data_stream);
    int get_no_clashes();

private:
    // Each of the two separate transformations applied to each unit.
    double alpha_z_left, beta_y_left, gamma_x_left, translate_x_left, translate_y_left, translate_z_left;
    double alpha_z_right, beta_y_right, gamma_x_right, translate_x_right, translate_y_right, translate_z_right;
    // Each pair of PDBs is transformed and then aligned as a whole on top of the calpha trace to get the RMSD.
    double rmsd;
    // Physics score
    soroban_score physics_score;
    // Number of clashed atoms (distance <3A)
    int no_clashes;
    string label_left, label_right;
};

#endif
