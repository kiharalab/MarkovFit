#include "pairwise_scores.h"
#include <fstream>
#include <utility>
//#include "align.h"
#include "atom.h"
#include "pairwise_features.h"
#include "point_transformation.h"
#include "rmsd.h"
#include "scoring.h"
// Included only because of the soroban_score predefined weights
#include "score_weights_options.h"
#include "soroban_score.h"
#include <sys/stat.h>

#ifdef WITH_MPI
// only include if the mpi version is being compiled
#include "mpi.h"
#endif

pairwise_scores::pairwise_scores(vector<string> pdb_unit_files, string calpha_trace_file, vector<string> unit_labels,
                         vector<pair<string, string> > neighbors, int num_transforms, string output_prefix, string pdb_output_prefix) :
unit_labels(unit_labels),
neighbors(neighbors),
num_transforms(num_transforms),
output_prefix(output_prefix), pdb_output_prefix(pdb_output_prefix) {
    // 1. Initialize the calpha trace, used to align and later to compute RMSDs.
    read_protein(calpha_trace_file, calpha_trace);
   
    // 2. Initialize Calpha fragments:
    pdb pdb_unit;
    pdb_unit.atoms.clear();
    for(size_t i = 0; i < calpha_trace.atoms.size(); i++) {
        if(((i+1) == calpha_trace.atoms.size())  || (calpha_trace.atoms[i].chain.compare(calpha_trace.atoms[i+1].chain)!=0 ) ) {
            pdb_unit.atoms.push_back(calpha_trace.atoms[i]);
            aligned_calphas.push_back(pdb_unit);
            pdb_unit.atoms.clear();
        } else {
            pdb_unit.atoms.push_back(calpha_trace.atoms[i]);
        }
    }
    
    // 3. Set the initial position for the PDB files (initial optimal placement).
    initialize_pdb_units(pdb_unit_files, unit_labels);
    // 4. For later use, initialize the mapping between labels and the unit index they represent
    for (size_t unit_index = 0; unit_index < unit_labels.size(); unit_index++) {
        string label = unit_labels[unit_index];
        label_indices[label] = unit_index;
    }
}

void pairwise_scores::initialize_pdb_units(vector<string> pdb_filenames, vector<string> labels) {
    for(size_t unit_index = 0; unit_index < pdb_filenames.size(); unit_index++) {
        pdb pdb_unit;
        read_protein(pdb_filenames[unit_index], pdb_unit);
        // Create the optimally aligned representation
        vector<atom> aligned_output;
        vector<atom> aligned_calpha;
        transformable_pdb aligned_unit(labels[unit_index], pdb_unit.atoms);
        aligned_unit.precompute_centroid();
        double centroid_x, centroid_y, centroid_z;
        aligned_unit.get_centroid(&centroid_x, &centroid_y, &centroid_z);
        aligned_units.push_back(aligned_unit);
    }
}

void pairwise_scores::generate_singleton_pdb(size_t unit_index, point_transformation_sequence& applied_transform,
                                         string pdb_prefix, size_t file_sequence_number) {
    transformable_pdb current_unit = aligned_units[unit_index];
    current_unit.apply_point_transformation_sequence(&applied_transform);
    write_complex(pdb_prefix + "-" + unit_labels[unit_index], file_sequence_number, current_unit.atoms);
}

void pairwise_scores::generate_feature_files(vector<string> transform_files) {
    // Ignore two body files if PDB's are expected. This is just to save time if only the PDB's are needed.
    if (pdb_output_prefix.empty()) {
        generate_two_body_files(transform_files, false);
    }
}

void pairwise_scores::generate_two_body_files(vector<string> transform_files, bool append_to_file) {
    vector <point_transformation> left_transforms, right_transforms;
    for (size_t edge_index = 0; edge_index < neighbors.size(); edge_index++) {
        pair<string,string> edge = neighbors[edge_index];
        size_t left_index = label_indices[edge.first];
        size_t right_index = label_indices[edge.second];
        left_transforms.clear();
        read_transformations(left_transforms, transform_files[left_index], edge.first);
        
        long long total_computations = left_transforms.size();
        
        if (total_computations < 0) {
            cerr << "Overflow in comparisons counter, density_lattice.get_significant_correlations" << endl;
        }
        string out;
#ifdef WITH_MPI
        // If mpi is used the different bases will be distributed between different ranks
        int int_rank, int_numtasks;
        MPI_Comm_rank(MPI_COMM_WORLD, &int_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &int_numtasks);
        // Have long long variables to force 64 bits. The number of total_computations
        // can be really large and these two variables are used in related calculations.
        long long rank = static_cast<long long>(int_rank);
        long long numtasks = static_cast<long long>(int_numtasks);
        long long iterations_per_task = total_computations / numtasks;
        long long total_remainder = total_computations % numtasks;
        long long remainder_processed = total_remainder > rank ? rank : total_remainder;
        long long start_index = iterations_per_task * rank + remainder_processed;
        long long end_index = start_index + iterations_per_task;
        // the remainder will be given 1 by 1 to each of the first ranks, until we run out of remaining bases
        if(rank < total_remainder)
        {
            end_index++;
        }
#else
        ios_base::openmode mode = append_to_file ? std::ofstream::app : std::ofstream::trunc;
        if (output_prefix.empty()) {
           out =  (unit_labels[left_index] + "-" + unit_labels[right_index] + ".mrf").c_str();
        } else {
            out =  (output_prefix + "-" + unit_labels[left_index] + "-" + unit_labels[right_index] + ".mrf").c_str();
        }
        ofstream sampling_results_file(out, mode);

        if (!append_to_file) {
            mrf_two_body_features::write_header(sampling_results_file);
        }
#endif
        right_transforms.clear();
        read_transformations(right_transforms, transform_files[right_index], edge.second);
        
        // Keep track of the iteration number with this variable
        long long current_iteration = 0;
        vector <mrf_two_body_features> feature_all;
        for(size_t i = 0; i<left_transforms.size(); i++) {
            current_iteration++;
#ifdef WITH_MPI
            if ((current_iteration - 1) < start_index || (current_iteration - 1) >= end_index) {
                // Ignore it, this rank is not supposed to compute it
                continue;
            }
#endif
            point_transformation base_left_transform = left_transforms[i];
            pair<size_t, point_transformation_sequence> left_samples = get_all_adjusted_transformations(left_index, base_left_transform);
            point_transformation_sequence applied_left_transform = left_samples.second;
            
            for(size_t j = 0; j<right_transforms.size(); j++) {
                point_transformation base_right_transform = right_transforms[j];
                pair<size_t, point_transformation_sequence> right_samples = get_all_adjusted_transformations(right_index,base_right_transform);
                point_transformation_sequence applied_right_transform = right_samples.second;
                mrf_two_body_features features = calculate_two_body_features(left_index, right_index, base_left_transform, base_right_transform,
                                                                             applied_left_transform, applied_right_transform);
#ifdef WITH_MPI
                feature_all.push_back(features);
#else
                features.write(sampling_results_file);
#endif
                
            }
        }
        //swap right_transforms
        vector<point_transformation>().swap( right_transforms );
        //swap left_transforms
        vector<point_transformation>().swap( left_transforms );
        
#ifdef WITH_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        int msg = rank;
        if (rank==0) {
            //write
            ios_base::openmode mode = append_to_file ? std::ofstream::app : std::ofstream::out;
            if (output_prefix.empty()) {
                out = (unit_labels[left_index] + "-" + unit_labels[right_index] + ".mrf").c_str();
            }else {
                out = (output_prefix + "-" + unit_labels[left_index] + "-" + unit_labels[right_index] + ".mrf").c_str();
            }
            ofstream sampling_results_file(out , mode);

            if (!append_to_file) {
                mrf_two_body_features::write_header(sampling_results_file);
            }
            for (size_t i=0; i<feature_all.size(); i++) {
                feature_all[i].write(sampling_results_file);
            }
            sampling_results_file.close();
            
            //send msg
            MPI_Send(&msg,1, MPI_INT, rank+1, rank, MPI_COMM_WORLD);
        } else {
            MPI_Recv(&msg,1,MPI_INT,rank-1, rank-1, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if (output_prefix.empty()) {
                out = (unit_labels[left_index] + "-" + unit_labels[right_index] + ".mrf").c_str();
            }else {
                out = (output_prefix + "-" + unit_labels[left_index] + "-" + unit_labels[right_index] + ".mrf").c_str();
            }

            ofstream sampling_results_file( out, std::ofstream::app);

            for (size_t i=0; i<feature_all.size(); i++) {
                feature_all[i].write(sampling_results_file);
            }
            sampling_results_file.close();
            if (rank+1<numtasks) {
                MPI_Send(&msg,1, MPI_INT, rank+1, rank, MPI_COMM_WORLD);
            }
        }
#else
        sampling_results_file.close();
#endif
    }
}

mrf_two_body_features pairwise_scores::calculate_two_body_features(size_t left_index, size_t right_index,
                                                               point_transformation& base_left_transform,
                                                               point_transformation& base_right_transform,
                                                               point_transformation_sequence& applied_left_transform,
                                                               point_transformation_sequence& applied_right_transform) {
    
    transformable_pdb left_pdb = aligned_units[left_index];
    transformable_pdb right_pdb = aligned_units[right_index];
    // Apply translations wrt the centroids
    left_pdb.apply_point_transformation_sequence(&applied_left_transform);
    right_pdb.apply_point_transformation_sequence(&applied_right_transform);
    pdb both_pdb;
    both_pdb.merge(left_pdb, right_pdb, both_pdb);

    vector<atom> input_matches;
    vector<vector<atom> > result;

    result =  both_pdb.get_all_matching_docked_calphas(aligned_calphas, &input_matches);

    double min =  DBL_MAX, rmsd;
    
    for(size_t j=0 ; j<result.size(); j++) {
        rmsd = calculate_allatom_rmsd(input_matches, result[j]);
        if(rmsd < min) {
            min = rmsd;
        }
    }
    
    //swap result
    vector<vector<atom> >().swap( result );
    //swap input_matches
    vector<atom>().swap( input_matches );
    
    // Calculate the physics score terms (it requires two separate vectors)
    vector<vector<atom> > decoy;
    decoy.push_back(left_pdb.atoms);
    decoy.push_back(right_pdb.atoms);
    
    int no_clashes = get_no_clashing_atoms(decoy);
    soroban_score score;
    score = compute_energy(decoy, TemplateWeights[0]);
    //swap decoy
    vector<vector<atom> >().swap( decoy );
    
    mrf_two_body_features features(base_left_transform.get_alpha(), base_left_transform.get_beta(),
                                   base_left_transform.get_gamma(), base_left_transform.get_tx(),
                                   base_left_transform.get_ty(), base_left_transform.get_tz(),
                                   base_right_transform.get_alpha(), base_right_transform.get_beta(),
                                   base_right_transform.get_gamma(), base_right_transform.get_tx(),
                                   base_right_transform.get_ty(), base_right_transform.get_tz(), min,
                                   score, no_clashes, unit_labels[left_index], unit_labels[right_index]);
    return features;
}

pair<size_t, point_transformation_sequence> pairwise_scores::get_all_adjusted_transformations(size_t unit_index, point_transformation base_transform) {
    pair<size_t, point_transformation_sequence> all;
    point_transformation_sequence adjusted_sequence;
    add_centroid_adjusted_transformation(unit_index, base_transform, &adjusted_sequence);
    all = std::make_pair(unit_index, adjusted_sequence);
    return all;
}

void pairwise_scores::add_centroid_adjusted_transformation(size_t unit_index, point_transformation transform, point_transformation_sequence* sequence) {
    double centroid_x, centroid_y, centroid_z;
    aligned_units[unit_index].precompute_centroid();
    aligned_units[unit_index].get_centroid(&centroid_x, &centroid_y, &centroid_z);
    sequence->add_translation(-centroid_x, -centroid_y, -centroid_z);
    sequence->add(transform);
}

void pairwise_scores::read_transformations(vector<point_transformation> &transformations, string transform_file, string unit) {
    size_t unit_index = label_indices[unit];
    ifstream current_file;
    int count = 0;
    //Read all input files and populate the 2D array
    
    //Read file
    current_file.open(transform_file.c_str());
    if (!current_file.is_open()) {
        cerr << "Could not open file." <<endl;
    }
    
    double a, b, c, d, e, f, g, h, w, z;
    vector<vector<double> > features;
    vector<double> feature;
   
    while (!current_file.eof()) {
        if(count>=num_transforms) {
            break;
        }
        while (current_file >> w >> a >> b >> c >> d >> e >> f >> g >> h >> z ) {
            point_transformation tmp(f, e, d, a, b, c);
            transformations.push_back(tmp);
            feature.clear();
            feature.push_back(g);
            feature.push_back(h);
            features.push_back(feature);
            count++;
            if(count>=num_transforms) {
                break;
            }
        }
    }
    
    current_file.close();
    current_file.clear();
    
#ifdef WITH_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank != 0) {
        //swap features, feature
        vector<vector<double> >().swap( features );
        return ;
    }
#endif
    //get unit
    string output_file;
    if (output_prefix.empty()) {
        output_file = unit+".mrf";
    } else {
        output_file = output_prefix+"_"+unit+".mrf";
    }
    
        ofstream result_file (output_file.c_str(),ios::out); //fix
        
        if (!result_file.is_open()) {
            cerr << "Could not create result file!" << endl;
            throw 1;
        }
        
        result_file << "RMSD\tCC\tOverlap\tAlpha_z\tBeta_y\tGamma_x\ttx\tty\ttz\tUnit"<< endl;

        //Write unit.mrf file
        for(size_t i=0; i< features.size(); i++) {
            transformable_pdb unit_pdb = aligned_units[unit_index];
            point_transformation base_transform = transformations[i];
            pair<size_t, point_transformation_sequence> samples = get_all_adjusted_transformations(unit_index,base_transform);
            point_transformation_sequence applied_transform = samples.second;
            
            unit_pdb.apply_point_transformation_sequence(&applied_transform);
            
            double centroid_x, centroid_y, centroid_z;
            unit_pdb.precompute_centroid();
            unit_pdb.get_centroid(&centroid_x, &centroid_y, &centroid_z);
            
            vector<atom> input_matches;
            vector<vector<atom> > result;
            result = unit_pdb.get_all_matching_calphas_multiple(aligned_calphas, &input_matches);
            
            double min =  DBL_MAX, rmsd;
            
            for(size_t j=0 ; j<result.size(); j++) {
                rmsd = calculate_allatom_rmsd(input_matches, result[j]);
                if(rmsd < min) {
                    min = rmsd;
                }
            }
            //swap input_matches
            vector<atom>().swap( input_matches );
            //swap result
            vector<vector<atom> >().swap( result );
            result_file << min << "\t" << features[i][0] << "\t" << features[i][1] << "\t" << transformations[i].get_alpha() << "\t"
            << transformations[i].get_beta() << "\t" << transformations[i].get_gamma() << "\t" << transformations[i].get_tx() << "\t"
            << transformations[i].get_ty() << "\t"<< transformations[i].get_tz() << "\t" << unit << endl;
        }
        result_file.close();
    vector<vector<double> >().swap( features );
}
