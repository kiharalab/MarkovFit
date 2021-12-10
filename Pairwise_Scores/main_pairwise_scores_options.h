#ifndef _MAIN_PDB_ALIGN_OPTIONS_H_
#define _MAIN_PDB_ALIGN_OPTIONS_H_

#include <getopt.h>
#include <cstdlib>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace std;

#define INPUT_PDBS_OPTION "input-pdb"
#define INPUT_PDBS_CHAR 'i'
#define INPUT_PDBS_DESCRIPTION "n full atom PDB files that are sampled around the centers of the c-alpha\n\t\ttemplate centroids"

#define TEMPLATE_CALPHA_PDB_OPTION "calpha"
#define TEMPLATE_CALPHA_PDB_CHAR 'c'
#define TEMPLATE_CALPHA_PDB_DESCRIPTION "c-alpha only PDB template file that represents the best fitting\n\t\tof the units in the complete map. It should contain c-alphas for all subunits in a single file."

#define LABELS_OPTION "labels"
#define LABELS_CHAR 'l'
#define LABELS_DESCRIPTION "Each of the n subunits are assigned a label. This comma-separated list of\n\t\ttext names is matched 1 to 1 in the same order as the full atom PDB files provided."

#define OUTPUT_PREFIX_OPTION "output-prefix"
#define OUTPUT_PREFIX_CHAR 'o'
#define OUTPUT_PREFIX_DESCRIPTION "Output filenames will start with this prefix"

#define PDB_OUTPUT_PREFIX_OPTION "pdb-output-prefix"
#define PDB_OUTPUT_PREFIX_CHAR 'p'
#define PDB_OUTPUT_PREFIX_DESCRIPTION "Triggers the output of PDB files that represent each singleton.\n\tOnly singleton sample files are generated if provided"

#define TRANSFORMATIONS_FILE_OPTION "transforms-file"
#define TRANSFORMATIONS_FILE_CHAR 't'
#define TRANSFORMATIONS_FILE_DESCRIPTION "Each of the n subunits have transformation file."

#define NUM_TRANSFORMS_OPTION "transforms-num"
#define NUM_TRANSFORMS_CHAR 'm'
#define NUM_TRANSFORMS_DESCRIPTION "Each of the n subunits have transformation file."

#define SHORT_OPTIONS ""

static struct option main_pairwise_scores_long_options[] =  {
    {INPUT_PDBS_OPTION, required_argument, 0, INPUT_PDBS_CHAR},
    {TEMPLATE_CALPHA_PDB_OPTION, required_argument, 0, TEMPLATE_CALPHA_PDB_CHAR},
    {LABELS_OPTION, required_argument, 0, LABELS_CHAR},
    {OUTPUT_PREFIX_OPTION, required_argument, 0, OUTPUT_PREFIX_CHAR},
    {PDB_OUTPUT_PREFIX_OPTION, required_argument, 0, PDB_OUTPUT_PREFIX_CHAR},
    {TRANSFORMATIONS_FILE_OPTION, required_argument, 0, TRANSFORMATIONS_FILE_CHAR},
    {NUM_TRANSFORMS_OPTION, required_argument, 0, NUM_TRANSFORMS_CHAR},
    {0,0,0,0}
};


class main_pairwise_scores_options {
public:
    main_pairwise_scores_options(int argc, char** argv) {
        labels_arg = template_pdb = output_prefix = pdb_output_prefix = "";
        num_transforms= 0;
        // so far we haven't encountered errors
        parse_error = false;
        bool options_left = true;
        
        while(options_left) {
            int option_index, c;
            c = getopt_long(argc, argv, SHORT_OPTIONS, main_pairwise_scores_long_options, &option_index);
            if(c == -1) {
                options_left = false;
            }
            else {
                // determine what type of option we got
                switch(c) {
                    case INPUT_PDBS_CHAR:
                        input_pdbs.push_back(optarg);
                        break;
                        
                    case TRANSFORMATIONS_FILE_CHAR:
                        transform_files.push_back(optarg);
                        break;
                        
                    case TEMPLATE_CALPHA_PDB_CHAR:
                        template_pdb = optarg;
                        break;
                        
                    case LABELS_CHAR:
                        labels_arg = optarg;
                        break;
                        
                    case NUM_TRANSFORMS_CHAR:
                        num_transforms = atof(optarg);
                        break;
                        
                    case OUTPUT_PREFIX_CHAR:
                        output_prefix = optarg;
                        break;
                        
                    case PDB_OUTPUT_PREFIX_CHAR:
                        pdb_output_prefix = optarg;
                        break;
                        
                    case '?': //option not recognized (an error is automatically printed)
                    default:
                        options_left = false;
                        parse_error = true;
                        break;
                }
            }
        }
    }
    
    // the parse was successful if the errors flag is false and if all required arguments were supplied
    bool parse_successful() {
        if(parse_error || labels_arg.empty() ||
           template_pdb.empty() || transform_files.empty()) {
            return false;
        }
        // Check if graph connections belong to label names
        vector<pair<string,string> > neighbors = get_neighbors();
        vector<string> labels = get_labels();
        set<string> labels_set(labels.begin(), labels.end());
        for(vector<pair<string,string> >::iterator it = neighbors.begin(); it != neighbors.end(); it++) {
            pair<string,string> neighbor_pair = *it;
            if(labels_set.find(neighbor_pair.first) == labels_set.end() ||
               labels_set.find(neighbor_pair.second) == labels_set.end()) {
                return false;
            }
        }
        // Check if the size of the input PDBs and EM maps are the same
        if((input_pdbs.size() != labels.size()) && (input_pdbs.size() != transform_files.size())) {
            return false;
        }
        return true;
    }
    
    string usage() {
        return string("Sample usage: pairwise_scores --input-pdb A.pdb --input-pdb B.pdb --input-pdb C.pdb ") +
        "-t ABC-calpha.pdb --transforms-file A_transform-file --transforms-file B_transform-file --transforms-file C_transform-file\n"
        "\t\t--labels A,B,C --calpha main_pdb --transforms-num 100 --output-prefix prefix \n\n"
        "Parameters:\n\t" +
        string(1, INPUT_PDBS_CHAR) + ": " + INPUT_PDBS_DESCRIPTION + "\n\t" +
        TRANSFORMATIONS_FILE_OPTION + ": " + TRANSFORMATIONS_FILE_DESCRIPTION + "\n\t" +
        LABELS_OPTION + ": " + LABELS_DESCRIPTION + "\n\t" +
        string(1, TEMPLATE_CALPHA_PDB_CHAR) + ": " + TEMPLATE_CALPHA_PDB_DESCRIPTION + "\n\t" +
        NUM_TRANSFORMS_OPTION + ": " + NUM_TRANSFORMS_DESCRIPTION + "\n\t" +
        OUTPUT_PREFIX_OPTION + ": " + OUTPUT_PREFIX_DESCRIPTION + "\n";
    }
    
    vector<string> get_input_pdbs() {
        return input_pdbs;
    }
    
    vector<string> get_transform_files() {
        return transform_files;
    }
    
    string get_template_pdb() {
        return template_pdb;
    }
    
    vector<string> get_labels() {
        vector<string> labels;
        string labels_param = labels_arg;
        while(!labels_param.empty())
        {
            size_t separator_pos = labels_param.find_first_of(',');
            string new_label = labels_param.substr(0, separator_pos);
            labels.push_back(new_label);
            labels_param = separator_pos ==  string::npos ? "" : labels_param.substr(separator_pos + 1);
        }
        return labels;
    }

    vector<pair<string,string> > get_neighbors() {
        vector<string> labels = get_labels();
        vector<pair<string, string> > neighbors;

        for(size_t i=0; i<labels.size()-1;i++) {
            for(size_t j=i+1; j<labels.size(); j++) {
                neighbors.push_back(make_pair(labels[i],labels[j]));
            }
        }

        return neighbors;
    }
    
    int get_num_transforms() {
        return num_transforms;
    }
    
    string get_output_prefix() {
        return output_prefix;
    }
    
    string get_pdb_output_prefix() {
        return pdb_output_prefix;
    }
    
private:
    vector<string> input_pdbs;
    vector<string> transform_files;
    int num_transforms;
    string template_pdb;
    string labels_arg;
    string output_prefix;
    string pdb_output_prefix;
    // set to false in case that unrecognized options or errors are found by getopt_long
    bool parse_error;
};

#endif
