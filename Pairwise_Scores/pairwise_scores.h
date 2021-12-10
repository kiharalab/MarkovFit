#ifndef _PAIRWISE_SCORES_H_
#define _PAIRWISE_SCORES_H_

#include <map>
#include <utility>
#include <vector>
#include "pairwise_features.h"
#include "pdb.h"
#include "transformable_pdb.h"

using std::pair;
using std::vector;

// Instances of this class are used to generate data files that sample different positions of PDB units fittted into an EM map.
// These files will contain features that are later used to train models using an MRF framework.
class pairwise_scores {
public:
    // pdb_unit_files: Each atomic description of the units e.g. A.pdb, B.pdb, etc
    // transform_files: Transformations to be applied to each subunit.
    // calpha_trace_file: Best arrangement of the PDB units that can be fitted into the complete_map_file.
    // All the pdb_unit_files can be superimposed onto this description of the c-alpha trace
    // unit_labels: Used to later identify the neighbor relations. Must be in the same order as the PDB and EM units provided
    // output_prefix: adds a prefix to sample files
    // pdb_output_prefix: if non-empty triggers the generation of PDB files (with filenames having this prefix) for each singleton transformation.
    //   If provided, pairwise files are not generated.
    pairwise_scores(vector<string> pdb_unit_files, string calpha_trace_file,
                vector<string> unit_labels, vector<pair<string, string> > neighbors, int num_transforms,
                string output_prefix, string pdb_output_prefix);
    // Main function that starts the process of generating the two-body feature files.
    void generate_feature_files(vector<string> transform_files);
private:
    void initialize_pdb_units(vector<string> pdb_filenames, vector<string> labels);
    // Generate a single PDB file by modifying protein[unit_index] with the transformation given and storing it in
    // <pdb_prefix>-<unit_label>-<file_sequence_number>.pdb
    void generate_singleton_pdb(size_t unit_index, point_transformation_sequence& applied_transformation,
                                string pdb_prefix, size_t file_sequence_number);
    // Computes the features for a particular transformation of the 2 vertices in the graph
    void generate_two_body_files(vector<string> transform_files, bool append_to_file);
    // correctly initialize the mrf_two_body_features instance.
    mrf_two_body_features calculate_two_body_features( size_t left_index, size_t right_index, point_transformation& base_left_transformation,
                                                      point_transformation& base_right_transformation,
                                                      point_transformation_sequence& applied_left_transformation,
                                                      point_transformation_sequence& applied_right_transformation);
    pair<size_t, point_transformation_sequence>  get_all_adjusted_transformations(size_t unit_index, point_transformation base_transformation);
    // Since transformations are calculated from each unit's centroid, this method adds three point_transformations to the sequence
    void add_centroid_adjusted_transformation(size_t unit_index, point_transformation transformation, point_transformation_sequence* sequence);
    void read_transformations(vector<point_transformation> &transformations, string transform_file, string unit);

    // It should be a PDB ID referenced from the EMDB entry that provides the best fitted superimposition of the c-alpha atoms into the map
    pdb calpha_trace;
    // Matches calpha_trace, loaded from an EM map obtained directly from the EMDB.
    // Full atomic PDB fragments that have been transformed to align as well as possible to the c-alpha trace.
    vector<transformable_pdb> aligned_units;
    vector<pdb> aligned_calphas;
    vector<point_transformation> alignment_transformations;
    vector<string> unit_labels;
    vector<string> transform_files;
    vector<pair<string,string> > neighbors;
    // Quick reference to determine what unit index corresponds to each label
    map<string, size_t> label_indices;
    int num_transforms;
    // Files created have this prefix.
    string output_prefix;
    // Same for PDB files
    string pdb_output_prefix;
};
#endif
