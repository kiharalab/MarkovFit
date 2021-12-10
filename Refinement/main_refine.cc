#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

#include "main_refine_options.h"
#include "refine.h"
#ifdef WITH_MPI
// only include if the mpi version is being compiled
#include "mpi.h"
#endif

int main(int argc, char** argv)
{
    main_refine_options options(argc,argv);
    if(!options.parse_successful())
    {
        std::cerr << options.usage();
        exit(EXIT_FAILURE);
    }
    
#ifdef WITH_MPI
    /*** Perform the MPI Initialization that applies to all ranks ***/
    MPI_Init(&argc,&argv);
#endif
    
    vector<string> input_pdb_files = options.get_input_pdbs();
    vector<string> transform_files = options.get_transform_files();
    string calpha_template_file = options.get_template_pdb();
    vector<string> labels = options.get_labels();
    vector<pair<string,string> > neighbors = options.get_neighbors();
    string output_prefix = options.get_output_prefix();
    string pdb_output_prefix = options.get_pdb_output_prefix();
    int num_transforms = options.get_num_transforms();
    
    refine sampler(input_pdb_files, calpha_template_file, labels,
                        neighbors, num_transforms, output_prefix, pdb_output_prefix);
        
    sampler.generate_feature_files(transform_files);

#ifdef WITH_MPI
    MPI_Finalize();
#endif

    exit(EXIT_SUCCESS);
}
