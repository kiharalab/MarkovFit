#ifndef _REFINE_H_
#define _REFINE_H_

#include <map>
#include <utility>
#include <vector>
#include "mrf_features.h"
#include "pdb.h"
#include "transformable_pdb.h"

using std::pair;
using std::vector;

// Instances of this class are used to generate data files that sample different positions of PDB units fittted into an EM map.
// These files will contain features that are later used to train models using an MRF framework.
class refine {
public:
    // pdb_unit_files: Each atomic description of the units e.g. A.pdb, B.pdb, etc
    // transform_files: Transformations to be applied to each subunit.
    // calpha_trace_file: Best arrangement of the PDB units that can be fitted into the complete_map_file.
    // output_prefix: adds a prefix to files
    refine(vector<string> pdb_unit_files, string calpha_trace_file,
                vector<string> unit_labels, vector<pair<string, string> > neighbors, int num_transforms,
                string output_prefix, string pdb_output_prefix);
    void read_transformations(vector<point_transformation> &transformations, string transform_file, string unit);

    // Main function that starts the process of generating the two-body feature files.
    void generate_feature_files(vector<string> transform_files);

private:
    void initialize_pdb_units(vector<string> pdb_filenames, vector<string> labels, pdb& calpha_template);
    void generate_singleton_pdb(size_t unit_index, point_transformation_sequence& applied_transformation,
                                string pdb_prefix, size_t file_sequence_number);
    // Computes the features for a particular transformation of the 2 subunits
    void generate_two_body_files(vector<string> transform_files, bool append_to_file);
    mrf_two_body_features calculate_two_body_features( size_t left_index, size_t right_index, point_transformation& base_left_transformation,
                                                      point_transformation& base_right_transformation,
                                                      point_transformation_sequence& applied_left_transformation,
                                                      point_transformation_sequence& applied_right_transformation);
    // Return value: a vector with one pair, where the unit_index is the same as the one input.
    pair<size_t, point_transformation_sequence>  get_all_adjusted_transformations(size_t unit_index, point_transformation base_transformation);
    // Since transformations are calculated from each unit's centroid, this method adds three point_transformations to the sequence:
    void add_centroid_adjusted_transformation(size_t unit_index, point_transformation transformation, point_transformation_sequence* sequence);
    
    void read_random_transformations(vector<point_transformation> &transformations, string transform_file);    
    vector<point_transformation> create_random_transformations(vector<point_transformation> transformations, string file_to_read, string file_to_write, string unit, int flag);

    void write_transformations(vector<point_transformation> transformations, string unit);
    void write_pairwise_transformations(vector<point_transformation> left_transformations,vector<point_transformation> right_transformations, string left_unit, string right_unit);
    void generate_one_body_files(double alpha_center,
                                              double beta_center, double gamma_center, double alpha_plus_minus,
                                              double beta_plus_minus, double gamma_plus_minus,
                                              double translation_plus_minus, double rotation_step,
                                 double translation_step, bool append_to_file, bool swapping_enabled);
    // It should be a PDB ID referenced from the EMDB entry that provides the best fitted superimposition of the c-alpha atoms into the map
    pdb calpha_trace;
    // Full atomic PDB fragments that have been transformed to align as well as possible to the c-alpha trace.
    vector<transformable_pdb> aligned_units;
    vector<pdb> aligned_calphas;
    // Transformations applied to input PDBs to optimally align them to the calpha_trace
    vector<point_transformation> alignment_transformations;
    vector<string> unit_labels;
    vector<string> transform_files;
    vector<pair<string,string> > neighbors;
    map<string, size_t> label_indices;
    int num_transforms;
    // Files created have this prefix.
    string output_prefix;
    // Same for PDB files
    string pdb_output_prefix;
};
#endif
