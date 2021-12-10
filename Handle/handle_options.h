#ifndef _HANDLE_OPTIONS_H_
#define _HANDLE_OPTIONS_H_

#include <cstdlib>
#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <cfloat>

using namespace std;

#define X_OPTION "correct-x"
#define X_CHAR 'x'
#define X_DESCRIPTION "X value of correct postion"

#define Y_OPTION "correct-y"
#define Y_CHAR 'y'
#define Y_DESCRIPTION "Y value of correct postion"

#define Z_OPTION "correct-z"
#define Z_CHAR 'z'
#define Z_DESCRIPTION "Z value of correct postion"

#define DISTANCE_OPTION "dist-threshold"
#define DISTANCE_CHAR 'd'
#define DISTANCE_DESCRIPTION "Distance threshold: < 5 Angstrom"

#define MIN_DISTANCE_OPTION "min-dist"
#define MIN_DISTANCE_CHAR 'm'
#define MIN_DISTANCE_DESCRIPTION "Minimum distance in Angstrom"

#define INPUT_OPTION "input-file"
#define INPUT_CHAR 'i'
#define INPUT_DESCRIPTION "Search result file"

#define OUTPUT_OPTION "output"
#define OUTPUT_CHAR 'o'
#define OUTPUT_DESCRIPTION "File where results will be stored"

// No short options in this binary
#define SHORT_OPTIONS ""

static struct option long_options[] = 
{
    {X_OPTION, required_argument, 0, X_CHAR},
    {Y_OPTION, required_argument, 0, Y_CHAR},
    {Z_OPTION, required_argument, 0, Z_CHAR},
    {DISTANCE_OPTION, required_argument, 0, DISTANCE_CHAR},
    {MIN_DISTANCE_OPTION, required_argument, 0, MIN_DISTANCE_CHAR},
    {INPUT_OPTION, required_argument, 0, INPUT_CHAR},
    {OUTPUT_OPTION, required_argument, 0, OUTPUT_CHAR},
    {0,0,0,0}
};

// This class parses the command line arguments supplied to the the 
// result handling program
class options
{
private:
    // These variables hold the parameters that come directly from the command line further processing is done by the methods
    //        int sort;
    float x,y,z,dist_threshold, min_dist;
    string output_filename;
    string input_filename;
    // set to false in case that unrecognized options or errors are found by getopt_long
    bool parse_error;
    
public:
    options(int argc, char** argv)
    {
        output_filename = "";
        input_filename = "";
        // default values, it's OK if they're not provided
        x = y = z = FLT_MAX;
        dist_threshold = 10;
        min_dist = 8;
        // so far we haven't encountered errors
        parse_error = false;
        
        bool options_left = true;
        
        while(options_left)
        {
            int option_index, c; // c represents the character associated with an option
            c = getopt_long(argc, argv, SHORT_OPTIONS, long_options, &option_index);
            if(c == -1)
            {
                options_left = false;
            }
            else
            {
                // determine what type of option we got
                switch(c)
                {
                    case X_CHAR:
                        x = atof(optarg);
                        break;
                        
                    case Y_CHAR:
                        y = atof(optarg);
                        break;
                        
                    case Z_CHAR:
                        z = atof(optarg);
                        break;
                        
                    case DISTANCE_CHAR:
                        dist_threshold = atof(optarg);
                        break;
                        
                    case MIN_DISTANCE_CHAR:
                        min_dist = atof(optarg);
                        break;
                        
                    case INPUT_CHAR:
                        input_filename = optarg;
                        break;
                        
                    case OUTPUT_CHAR:
                        output_filename = optarg;
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
    
    // the parse was successful if the errors flag is false and if all required
    // arguments were supplied
    bool parse_successful()
    {
        return !parse_error && !output_filename.empty() && !input_filename.empty();
    }
    
    string usage()
    {
        return
        string("Usage: handle  ") +
        "[--correct-x x_value] [--correct-y y_value] " +
        "[--correct-z z_value]  "  + "[--dist-threshold distance_threshold]  " +
        "[--min-dist min dist in Angstrom] "+
        "[--input-file input file name]  " + " --output output file name\n\n" +
        "Parameters:\n\t" +
        X_OPTION + ": " + X_DESCRIPTION + "\n\t" +
        Y_OPTION + ": " + Y_DESCRIPTION + "\n\t" +
        Z_OPTION + ": " + Z_DESCRIPTION + "\n\t" +
        DISTANCE_OPTION + ": " + DISTANCE_DESCRIPTION + "\n\t" +
        MIN_DISTANCE_OPTION + ": " + MIN_DISTANCE_DESCRIPTION + "\n\t" +
        INPUT_OPTION + ": " + INPUT_DESCRIPTION + "\n\t" +
        OUTPUT_OPTION + ": " + OUTPUT_DESCRIPTION + "\n";
    }
    
    // Value of x axes of correct position
    float get_x_value()
    {
        return x;
    }
    
    // Value of y axes of correct position
    float get_y_value()
    {
        return y;
    }
    
    // Value of z axes of correct position
    float get_z_value()
    {
        return z;
    }
    
    // Value of distance threshold
    float get_dist_threshold()
    {
        return dist_threshold;
    }
    
    // min distance for clustering
    float get_min_dist()
    {
        return min_dist;
    }

    // Path to the output file created
    string get_output_filename()
    {
        return output_filename;
    }
    
    // input result files to analyze
    string get_input_filename()
    {
        return input_filename;
    }
    
};

#endif
