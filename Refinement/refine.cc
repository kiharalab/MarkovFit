#include "refine.h"
#include <fstream>
#include <utility>
#include "atom.h"
#include "mrf_features.h"
#include "point_transformation.h"
#include "rmsd.h"
#include "scoring.h"
#include "score_weights_options.h"
#include "soroban_score.h"
#include <sys/stat.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>

#ifdef WITH_MPI
// only include if the mpi version is being compiled
#include "mpi.h"
#endif

refine::refine(vector<string> pdb_unit_files, string calpha_trace_file, vector<string> unit_labels,
                         vector<pair<string, string> > neighbors, int num_transforms, string output_prefix, string pdb_output_prefix) :
unit_labels(unit_labels),
neighbors(neighbors),
num_transforms(num_transforms),
output_prefix(output_prefix), pdb_output_prefix(pdb_output_prefix) {
    srand(time(0));
    // 1. Initialize the calpha trace, used to align and later to compute RMSDs.
    read_protein(calpha_trace_file, calpha_trace, true);
    
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
    initialize_pdb_units(pdb_unit_files, unit_labels, calpha_trace);
    // 4. For later use, initialize the mapping between labels and the unit index they represent
    for (size_t unit_index = 0; unit_index < unit_labels.size(); unit_index++) {
        string label = unit_labels[unit_index];
        label_indices[label] = unit_index;
    }
}

void refine::initialize_pdb_units(vector<string> pdb_filenames, vector<string> labels, pdb& calpha_template) {
    for(size_t unit_index = 0; unit_index < pdb_filenames.size(); unit_index++) {
        pdb pdb_unit;
        read_protein(pdb_filenames[unit_index], pdb_unit,false);
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

void refine::generate_singleton_pdb(size_t unit_index, point_transformation_sequence& applied_transform,
                                         string pdb_prefix, size_t file_sequence_number) {
    transformable_pdb current_unit = aligned_units[unit_index];
    current_unit.apply_point_transformation_sequence(&applied_transform);
    write_complex(pdb_prefix + "-" + unit_labels[unit_index], file_sequence_number, current_unit.atoms);
}

void refine::generate_feature_files(vector<string> transform_files) {
    // Ignore two body files if PDB's are expected. This is just to save time if only the PDB's are needed.
    if (pdb_output_prefix.empty()) {
        generate_two_body_files(transform_files, false);
    }
}
/////////////////////////////////////////////////////////////
void refine::generate_two_body_files(vector<string> transform_files, bool append_to_file) {
    vector <point_transformation> left_transforms, right_transforms;
    
    int int_rank, int_numtasks;
    MPI_Comm_rank(MPI_COMM_WORLD, &int_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &int_numtasks);
    string file_to_write = "random_transforms.txt";
    int no_transf = 10, no_random = 50;
    vector<vector<point_transformation> > all_transf;
    all_transf.clear();
    
    //Process 0 reads transformations from files
    if (int_rank==0) {
        //create transformations
        for (size_t i=0; i<unit_labels.size(); i++) {
            string chain_id = unit_labels[i];
            string transform_file = chain_id+"_"+file_to_write;
            vector <point_transformation> chain_transforms;
            chain_transforms.clear();
            read_random_transformations(chain_transforms, transform_files[i]);
            all_transf.push_back(chain_transforms);
            vector <point_transformation> transformations;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    for (int rot = 0; rot < 10; rot++) {
        //Process 0 generate random transformations and sends them to other processes
        cout<< "Rotation # "<< rot+1 <<" :: Process "<<int_rank << endl;
        if (int_rank==0) {
            //create transformations
            for (size_t i=0; i<unit_labels.size(); i++) {
                string chain_id = unit_labels[i];
                string transform_file = chain_id + "_" + file_to_write;
                
                vector <point_transformation> chain_transforms;
                chain_transforms.clear();
                chain_transforms = all_transf[i];
                create_random_transformations(chain_transforms, transform_file, transform_file, chain_id, 2);
            }
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        all_transf.clear();
        //Read transformations
        for (size_t i=0; i<unit_labels.size(); i++) {
            string chain_id = unit_labels[i];
            string transform_file = chain_id+"_"+file_to_write;
            vector <point_transformation> chain_transforms;
            chain_transforms.clear();
            read_random_transformations(chain_transforms, transform_file);
            all_transf.push_back(chain_transforms);
        }
        //Compute physics score
        long long total_computations = no_transf;
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
        
        // Keep track of the iteration number with this variable
        vector <mrf_two_body_features> feature_all;
        
        //Loop over all transformations
        for (size_t i = 0; i < unit_labels.size(); i++) {
            double sum[no_random*no_transf] = {0};

            for (int j = 0; j < no_transf; j++) {
#ifdef WITH_MPI
                if (j < start_index || j >= end_index) {continue;}
#endif
                for (size_t edge_index = 0; edge_index < neighbors.size(); edge_index++) {
                    pair<string,string> edge = neighbors[edge_index];
                    size_t left_index = label_indices[edge.first];
                    size_t right_index = label_indices[edge.second];
                    if (left_index != i && right_index != i) {continue;}

                    for (int k=0; k<no_random; k++) {
                        if (left_index == i) {
                            point_transformation base_left_transform = all_transf[left_index][j*no_random+k];
                            pair<size_t, point_transformation_sequence> left_samples =
                            get_all_adjusted_transformations(left_index, base_left_transform);
                            point_transformation_sequence applied_left_transform = left_samples.second;
                            
                            point_transformation base_right_transform = all_transf[right_index][j*no_random];
                            pair<size_t, point_transformation_sequence> right_samples = get_all_adjusted_transformations(right_index,base_right_transform);
                            point_transformation_sequence applied_right_transform = right_samples.second;
                            mrf_two_body_features features = calculate_two_body_features(left_index, right_index,
                                                                                         base_left_transform, base_right_transform, applied_left_transform, applied_right_transform);
                            
                            sum[j*no_random+k] += (-1*features.get_physics_score());

                        }

                        if (right_index == i) {
                            point_transformation base_left_transform = all_transf[left_index][j];
                            pair<size_t, point_transformation_sequence> left_samples =
                            get_all_adjusted_transformations(left_index, base_left_transform);
                            point_transformation_sequence applied_left_transform = left_samples.second;
                            
                            point_transformation base_right_transform = all_transf[right_index][j*no_random+k];
                            pair<size_t, point_transformation_sequence> right_samples = get_all_adjusted_transformations(right_index,base_right_transform);
                            point_transformation_sequence applied_right_transform = right_samples.second;
                            mrf_two_body_features features = calculate_two_body_features(left_index, right_index,
                                                                                         base_left_transform, base_right_transform, applied_left_transform, applied_right_transform);
                            
                            sum[j*no_random+k] += (-1*features.get_physics_score());
                        }
                    }//End random
                }//End Edges
            }//End Transformations
            
            //Send sum to all other processes
            if (rank == 0) {
                for (int n=1; n<numtasks; n++) {
                    double temp_sum[no_random*no_transf]={0};
                    MPI_Recv(temp_sum,no_random*no_transf,MPI_DOUBLE,n, n, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                    for (int k=0; k<no_random*no_transf; k++) {
                        sum[k]+=temp_sum[k];
                    }
                }
                
            } else {
                MPI_Send(sum,no_random*no_transf,MPI_DOUBLE,0, rank, MPI_COMM_WORLD);
            }
            
            //Process 0 sends final sum to other processes
            if (rank == 0) {
                for (int i=1; i<numtasks; i++) {
                    MPI_Send(sum,no_random*no_transf,MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                }
            } else {
                double temp_sum[no_random*no_transf]={0};
                MPI_Recv(temp_sum,no_random*no_transf,MPI_DOUBLE,0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);

                for (int k=0; k<no_random*no_transf; k++) {
                    sum[k]=temp_sum[k];
                }
            }
            
            double max[no_transf];
            int max_index[no_transf];
            
            for (int i=0; i<no_transf; i++) {
                max[i]=-99999999999999;
                max_index[i]=-1;
            }
            //Find best transformation
            for (int x=0; x<no_transf; x++) {
                for (int j=0; j<no_random; j++) {
                    if (sum[x*no_random+j]>max[x]) {
                        max[x] = sum[x*no_random+j];
                        max_index[x]= x*no_random+j ;
                    }
                }
            }
            
            //Keep best transformations and remove all others
            int count = no_transf-1;
            for (int k=no_transf*no_random-1; k>=0; k--) {
                if (k!=max_index[count]) {
                    all_transf[i].erase(all_transf[i].begin()+k);
                } else {
                    count--;
                }
            }
        }//End Chain
    }
    if (int_rank==0) {
        for (size_t i=0; i<all_transf.size(); i++) {
            //cout<<"Chain "<<i<<" ("<<unit_labels[i]<<") \n";
            write_transformations(all_transf[i], unit_labels[i]);
            for (size_t j=i+1; j<all_transf.size(); j++) {
                write_pairwise_transformations(all_transf[i],all_transf[j], unit_labels[i],unit_labels[j]);
            }
        }
    }
    
}

mrf_two_body_features refine::calculate_two_body_features(size_t left_index, size_t right_index,
                                                               point_transformation& base_left_transform,
                                                               point_transformation& base_right_transform,
                                                               point_transformation_sequence& applied_left_transform,
                                                               point_transformation_sequence& applied_right_transform) {
    
    transformable_pdb left_pdb = aligned_units[left_index];
    transformable_pdb right_pdb = aligned_units[right_index];
    // Apply translations wrt the centroids
    left_pdb.apply_point_transformation_sequence(&applied_left_transform);
    right_pdb.apply_point_transformation_sequence(&applied_right_transform);
    
    double min = 0;
    // Calculate the physics score terms (it requires two separate vectors)
    vector<vector<atom> > decoy;
    decoy.push_back(left_pdb.atoms);
    decoy.push_back(right_pdb.atoms);
    int no_clashes = get_no_clashing_atoms(decoy);

    soroban_score score = compute_energy(decoy, TemplateWeights[0]);

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

pair<size_t, point_transformation_sequence> refine::get_all_adjusted_transformations(size_t unit_index, point_transformation base_transform) {
    pair<size_t, point_transformation_sequence> all;
    point_transformation_sequence adjusted_sequence;
    add_centroid_adjusted_transformation(unit_index, base_transform, &adjusted_sequence);
    all = std::make_pair(unit_index, adjusted_sequence);
    return all;
}

void refine::add_centroid_adjusted_transformation(size_t unit_index, point_transformation transform, point_transformation_sequence* sequence) {
    double centroid_x, centroid_y, centroid_z;
    aligned_units[unit_index].get_centroid(&centroid_x, &centroid_y, &centroid_z);
    sequence->add_translation(-centroid_x, -centroid_y, -centroid_z);
    sequence->add(transform);
}

void refine::read_random_transformations(vector<point_transformation> &transformations, string transform_file) {
    ifstream current_file;
    
    //Read file
    current_file.open(transform_file.c_str());
    if (!current_file.is_open()) {
        cerr << "Could not open file.." <<endl;
    } else {
    }
    
    double a, b, c, d, e, f, h;
    while (!current_file.eof()) {
        while (current_file >> h >>  a >> b >> c >> d >> e >> f ) {
            point_transformation tmp(a, b, c, d, e, f);
            transformations.push_back(tmp);
        }
    }
    
    current_file.close();
    current_file.clear();
}

vector<point_transformation> refine::create_random_transformations(vector<point_transformation> transformations, string file_to_read, string file_to_write, string unit, int flag) {
    
    if (flag == 1) {
        ifstream current_file;
        //Read file
        current_file.open(file_to_read.c_str());
        if (!current_file.is_open()) {
            cerr << "Could not open file..." <<endl;
        }
        double a, b, c, d, e, f, n;
        while (!current_file.eof()) {
            while (current_file >> n>> a >> b >> c >> d >> e >> f  ) {
                point_transformation tmp(d, e, f, a, b, c);
                transformations.push_back(tmp);
            }
        }
        current_file.close();
        current_file.clear();
    }
    
#ifdef WITH_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
    int upper, lower,num;
    vector<point_transformation> final_transformations;
    
    //create random transformations
    for (size_t i=0; i<transformations.size(); i++) {
        int num1, num2, num3;
        for (size_t j=0; j<50; j++) {
            point_transformation tmp;
            tmp = transformations[i];
            
            if (j==0) {
                tmp = transformations[i];
            } else if(j>=1 && j<=5) {
                upper = 2; lower = -2;
                num = (rand() %(upper - lower + 1)) + lower;
                tmp.set_translation(tmp.get_tx()+num,tmp.get_ty(),tmp.get_tz());
            } else if(j>=6 && j<=10) {
                upper = 2; lower = -2;
                num = (rand() %(upper - lower + 1)) + lower;
                tmp.set_translation(tmp.get_tx(),tmp.get_ty()+num,tmp.get_tz());
            } else if(j>=11 && j<=15) {
                upper = 2; lower = -2;
                num = (rand() %(upper - lower + 1)) + lower;
                tmp.set_translation(tmp.get_tx(),tmp.get_ty(),tmp.get_tz()+num);
            } else if(j>=16 && j<=20) {
                upper = 10; lower = -10;
                num = (rand() %(upper - lower + 1)) + lower;
                num += tmp.get_alpha();
                if (num > 360) { num -= 360; }
                tmp.set_rotation(num,tmp.get_beta(),tmp.get_gamma()); // rand(0,5)
            } else if(j>=21 && j<=25) {
                upper = 10; lower = -10;
                num = (rand() %(upper - lower + 1)) + lower;
                num += tmp.get_beta();
                if (num > 360) { num -= 360; }
                tmp.set_rotation(tmp.get_alpha(),num,tmp.get_gamma()); // rand(0,5)
            } else if(j>=26 && j<=30) {
                upper = 10; lower = -10;
                num = (rand() %(upper - lower + 1)) + lower;
                num += tmp.get_gamma();
                if (num > 360) { num -= 360; }
                tmp.set_rotation(tmp.get_alpha(),tmp.get_beta(),num); // rand(0,5)
            } else if(j>=30 && j<=35) {
                upper = 2; lower = -2;
                num1 = (rand() %(upper - lower + 1)) + lower;
                num2 = (rand() %(upper - lower + 1)) + lower;
                num3 = (rand() %(upper - lower + 1)) + lower;
                tmp.set_translation(tmp.get_tx()+num1,tmp.get_ty()+num2,tmp.get_tz()+num3); // rand(0,5)
            } else if(j>=36 && j<=40) {
                upper = 10; lower = -10;
                num1 = (rand() %(upper - lower + 1)) + lower;
                num2 = (rand() %(upper - lower + 1)) + lower;
                num3 = (rand() %(upper - lower + 1)) + lower;
                num1 += tmp.get_alpha();
                num2 += tmp.get_beta();
                num3 += tmp.get_gamma();
                if (num1 > 360) { num1 -= 360; }
                if (num2 > 360) { num2 -= 360; }
                if (num3 > 360) { num3 -= 360; }
                tmp.set_rotation(num1, num2, num3); // rand(0,5)
            } else {
                upper = 10; lower = -10;
                num1 = (rand() %(upper - lower + 1)) + lower;
                num2 = (rand() %(upper - lower + 1)) + lower;
                num3 = (rand() %(upper - lower + 1)) + lower;
                num1 += tmp.get_alpha();
                num2 += tmp.get_beta();
                num3 += tmp.get_gamma();
                if (num1 > 360) { num1 -= 360; }
                if (num2 > 360) { num2 -= 360; }
                if (num3 > 360) { num3 -= 360; }
                tmp.set_rotation(num1, num2, num3);
                
                upper = 2; lower = -2;
                num1 = (rand() %(upper - lower + 1)) + lower;
                num2 = (rand() %(upper - lower + 1)) + lower;
                num3 = (rand() %(upper - lower + 1)) + lower;
                
                tmp.set_translation(tmp.get_tx()+num1,tmp.get_ty()+num2,tmp.get_tz()+num3);
            }
            final_transformations.push_back(tmp);
        }
    }
    
    //get unit
    ofstream result_file (file_to_write.c_str()); //fix
    
    if (!result_file.is_open()) {
        cerr << "Could not create result file!" << endl;
        throw 1;
    }
    
    //Write file
    for(size_t i=0; i< final_transformations.size(); i++) {
        result_file << "0\t"<< final_transformations[i].get_alpha() << "\t" << final_transformations[i].get_beta()
        << "\t" << final_transformations[i].get_gamma() << "\t" << final_transformations[i].get_tx() << "\t" << final_transformations[i].get_ty() << "\t" << final_transformations[i].get_tz() << endl;
    }
    result_file.close();
    return final_transformations;
#endif
}

//Write Transformations
void refine::write_transformations(vector<point_transformation> transformations, string unit) {
    size_t unit_index = label_indices[unit];
    
    //get unit
    string output_file;
    if (output_prefix.empty()) {
        output_file = unit+".refine";
    } else {
        output_file = output_prefix+"_"+unit+".refine";
    }
    ofstream result_file (output_file.c_str());
    
    
    if (!result_file.is_open()) {
        cerr << "Could not create result file!" << endl;
        throw 1;
    }
    
    result_file << std::fixed << "RMSD\tAlpha_z\tBeta_y\tGamma_x\ttx\tty\ttz\tUnit"<< endl;
    
    //Write unit.mrf file
    for(size_t i=0; i< transformations.size(); i++) {
        transformable_pdb unit_pdb = aligned_units[unit_index];
        point_transformation base_transform = transformations[i];
        pair<size_t, point_transformation_sequence> samples = get_all_adjusted_transformations(unit_index,base_transform);
        point_transformation_sequence applied_transform = samples.second;
        
        unit_pdb.apply_point_transformation_sequence(&applied_transform);
        
        double centroid_x, centroid_y, centroid_z;
        unit_pdb.precompute_centroid();
        unit_pdb.get_centroid(&centroid_x, &centroid_y, &centroid_z);
        
        //Write pdb structure
        if (output_prefix.empty()) {
            write_complex(unit+"_phy_refine_decoy", i, unit_pdb.atoms);
        } else {
            write_complex(output_prefix+"_"+unit+"_phy_refine_decoy", i, unit_pdb.atoms);
        }
        
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
        result_file << std::fixed << min << "\t"<<  transformations[i].get_alpha() << "\t" << transformations[i].get_beta()
        << "\t" << transformations[i].get_gamma() << "\t" << transformations[i].get_tx() << "\t" << transformations[i].get_ty() << "\t"
        << transformations[i].get_tz() << "\t"  << unit << endl;
        
    }
    result_file.close();
}

//Write Pairwise Transformations
void refine::write_pairwise_transformations(vector<point_transformation> left_transformations,vector<point_transformation> right_transformations, string left_unit, string right_unit) {
    
    //get unit
    string output_file;
    if (output_prefix.empty()) {
        output_file = left_unit+"-"+right_unit+".refine";
    } else {
        output_file = output_prefix+"_"+left_unit+"-"+right_unit+".refine";
    }
    ofstream result_file (output_file.c_str());
    
    if (!result_file.is_open()) {
        cerr << "Could not create result file!" << endl;
        throw 1;
    }
    
    mrf_two_body_features::write_header(result_file);
    
    size_t left_index = label_indices[left_unit];
    size_t right_index = label_indices[right_unit];
    
    for (size_t i=0; i<left_transformations.size(); i++) {
        point_transformation base_left_transform = left_transformations[i];
        pair<size_t, point_transformation_sequence> left_samples =
        get_all_adjusted_transformations(left_index, base_left_transform);
        point_transformation_sequence applied_left_transform = left_samples.second;
        
        point_transformation base_right_transform = right_transformations[i];
        pair<size_t, point_transformation_sequence> right_samples = get_all_adjusted_transformations(right_index,base_right_transform);
        point_transformation_sequence applied_right_transform = right_samples.second;
        
        mrf_two_body_features features = calculate_two_body_features(left_index, right_index,
                                                                     base_left_transform, base_right_transform, applied_left_transform, applied_right_transform);
        features.write(result_file);
    }//End transformations
    result_file.close();
}

//void refine::read_transformations(vector<point_transformation> &transformations, string transform_file, string unit) {
//    ifstream current_file;
//    int count = 0;
//    //Read all input files and populate the 2D array
//    //Read file
//    current_file.open(transform_file.c_str());
//    if (!current_file.is_open()) {
//        cerr << "Could not open file." <<endl;
//    } else {
//        //cout << "Open file " << transform_file.c_str() <<"\n";
//    }
//    double a, b, c, d, e, f, g, h, w, z;
//    //int l, m ,n;
//    vector<vector<double> > features;
//    vector<double> feature;
//
//    while (!current_file.eof()) {
//        if(count>=num_transforms) {
//            break;
//        }
//        while (current_file >> w >> a >> b >> c >> d >> e >> f >> g >> h >> z ) {
//            point_transformation tmp(f, e, d, a, b, c);
//            transformations.push_back(tmp);
//            feature.clear();
//            feature.push_back(g);
//            feature.push_back(h);
//            features.push_back(feature);
//            count++;
//            if(count>=num_transforms) {
//                break;
//            }
//        }
//    }
//
//    current_file.close();
//    current_file.clear();
//
//#ifdef WITH_MPI
//    int rank;
//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//    if (rank != 0) {
//        return ;
//    }
//#endif
//}
