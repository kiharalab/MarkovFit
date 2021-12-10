#ifndef _SCORE_WEIGHTS_OPTIONS_H_
#define _SCORE_WEIGHTS_OPTIONS_H_

// This header files defines the string values used to parse command line
// options related to the weights applied to the 9-12 terms used for scoring

#include <fstream>
#include "utils.h"

// Constants used passed to getopt_long to help in the parsing process
#define WEIGHTS_OPTION "weights"
#define WEIGHTS_CHAR 'w'

// Constants that represent each of the weight configurations (these are the values given to the program in the command line args)
#define WEIGHT_TYPE_ALL string("all")
#define WEIGHT_TYPE_NO_OVERALL_NOSASA string("novdw-noelect-nosasa")
#define WEIGHT_TYPE_NO_OVERALL_NOSASA_NOLZSOLV string("novdw-noelect-nosasa-nolzsolv")
#define WEIGHT_TYPE_NO_OVERALL_NO_LZSOLV string("novdw-noelect-nolzsolv")

// description of the command line parameter (for display purposes)
#define WEIGHT_PARAM_DESCRIPTION "[--weights <" + WEIGHT_TYPE_ALL + "/" + WEIGHT_TYPE_NO_OVERALL_NOSASA + "/" + \
				WEIGHT_TYPE_NO_OVERALL_NOSASA_NOLZSOLV + "/" + WEIGHT_TYPE_NO_OVERALL_NO_LZSOLV + "/weightfile.txt>]"

// simple macro that exhaustively compares the string provided against each of the 4 possibilities and returns
// the appropriate value from the soroban_score::WeightType enum if it returns -1 then it's not a valid weight type
#define GET_WEIGHT_TYPE(input) ((!WEIGHT_TYPE_ALL.compare(input)) ? 0 : (!WEIGHT_TYPE_NO_OVERALL_NOSASA.compare(input)) ? 1 :\
										(!WEIGHT_TYPE_NO_OVERALL_NOSASA_NOLZSOLV.compare(input)) ? 2 :\
										(!WEIGHT_TYPE_NO_OVERALL_NO_LZSOLV.compare(input)) ? 3 : -1)

// Holds all the predefined sets of weights used
static double TemplateWeights[5][12] =
{
	// Weights used when considering all 12 terms
	{53.074, 151.294, 0.15, 1.373, 13.235, -0.097, 15.93, -0.734, 0, 0, 27.635, 79.063},
	// Weights used when using 9 terms, excluding overall van der waals, electrostatics and solvation
	{0, 81.9602, 0.1916921, 0, 10.57919, -0.05739181, 12.37243, -0.5037819, 0, 2.853197, 0, 58.7541},
	// The next line holds the original9 term weights
	//{0, 86.169, 0.12231 ,0, 9.8761 ,1.3775, 8.2947, 0.7689 ,1.4997 ,-40.51, 0, 86.205},
	// Reweighted without solvation terms because they were believed to be adding noise
	{0, 126.0055622, 0.47318348, 0, 13.03847991, -0.09138538, 16.64044339, -0.59111153, 0, 0, 0, 64.91418416},
	// Using 9 terms but Global solvation instead of Gaussian
	{0, 73.7674467196506, 0.192224163776118, 0, 10.4186923684025, -0.0234724295622823, 11.9451180592874,
		-0.542866981812019, 0, 0, 27.7843927854426, 42.1649879934922},
    // Taken from protein-protein docking
    {0, 0.338, 0.08, 0, 0.025, 0.002, 0.025, 0.098, 0.441, 0.279, 0.344, 0.164}

};

class score_weights_options
{

	public:
		/*
		 * Initializes the default type of weight used
		 */
		score_weights_options()
		{
			weights = TemplateWeights[1]; // this is the 9 term version with no overall values and no sasa
		}

		/*
		 * Returns a double-array of size 12 that contains the weights used to score
		 * It simply accesses the pointer stored by this instance which can be either the default one
		 * or the one specified through a call to set_weights
		 */
		double* get_weights()
		{
			return weights;
		}
		
	protected:
		/*
		 * This method should be invoked if a command line param specifies a change
		 * in the default weights used. The string expected is either one of the 4 types of descriptions
		 * given to the template weights or, if it doesn't match any of those, it will be interpreted as a
		 * filename that contains a custom set of weights
		 */
		void set_weights(string weight_description)
		{
			int weight_type = GET_WEIGHT_TYPE(weight_description);
			if(weight_type != -1) //it's one of the default values
			{
				weights = TemplateWeights[weight_type];
			}
			else // custom weights
			{
				ifstream weights_stream(weight_description.c_str(), ifstream::in);
				if(weights_stream.good())
				{
					string allweights;
					weights_stream >> allweights;
					vector<string> separateweights = split(allweights, ',');
					weights = new double[separateweights.size()];
					for(size_t weightindex = 0; weightindex < separateweights.size(); weightindex++)
					{
						weights[weightindex] = atof(separateweights[weightindex].c_str());
					}

					weights_stream.close();
				}
				else
				{
					// if it's not a file then set as null pointer
					weights = 0;
				}
			}
		}
	private:
		/*
		 * Variable that points to the correct set of weights to be used by the main program
		 */
		double* weights;

};

#endif
