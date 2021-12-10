#ifndef _SOROBAN_SCORE_H_
#define _SOROBAN_SCORE_H_

#include <iomanip>
#include <iostream>

using namespace std;

typedef struct soroban_t
{
	// van der waals score part
	// compute overall vdw and vdw as attractive and repulsive
	double vdw, vdw_attr, vdw_rep;
	// compute overall electrostatics
	// compute electrostatics which is split into
	// short range attractive elec_sr_attr
	// short range repulsive elec_sr_rep
	// long range attractive elec_lr_attr
	// long range repulsive elec_lr_rep
	double elec, elec_sr_attr, elec_sr_rep, elec_lr_attr, elec_lr_rep;
	// hydrogen bonding potential + disulphide
	double hbp_ss;
	// compute gaussian based solvation energy over all atoms
	// compute solvation energy based on accessible surface area
	double solv, sasa;
	// desolvation based on atom contact potential    
	double acp;
} soroban_t;

class soroban_score
{
	private:
		// Contains all 12 values that compose the score
		soroban_t score;
		// This enum is created to assure that the weight order will always be the same
		// Each of the different weight definitions used by the scoring scheme should follow this order
		enum WeightIndex
		{
			vdw_weight = 0,
			vdw_attr_weight = 1,
			vdw_rep_weight = 2,
			elec_weight = 3,
			elec_sr_attr_weight = 4,
			elec_lr_attr_weight = 5,
			elec_sr_rep_weight = 6,
			elec_lr_rep_weight = 7,
			hbp_ss_weight = 8,
			solv_weight = 9,
			sasa_weight = 10,
			acp_weight = 11
		};
	
		// Stores a pointer to the array containing the weights used to express
		// the final score as a weighted linear combination
		// of all the terms
		double* weights;
	public:

		// default constructor sets the weights to a NULL pointer
		soroban_score()
		{
			weights = NULL;
			score.vdw = score.vdw_attr = score.vdw_rep = score.elec = score.elec_sr_attr = score.elec_sr_rep =
				score.elec_lr_attr = score.elec_lr_rep = score.hbp_ss = score.solv = score.sasa = score.acp = 0;
		}

		// The parameter supplied describes the weights that should be applied to each term when treating the
		// overall score as a linear combination of all the terms.
		soroban_score(double* weights)
		{
			this->weights = weights;
			score.vdw = score.vdw_attr = score.vdw_rep = score.elec = score.elec_sr_attr = score.elec_sr_rep =
				score.elec_lr_attr = score.elec_lr_rep = score.hbp_ss = score.solv = score.sasa = score.acp = 0;
		}

		/*
		 * Provides access to the basic struct that contains the
		 * 12 scoring values
		 */
		soroban_t get_base()
		{
			return score;
		}

		/*
		 * Overwrites the internal structure that holds the 12 scoring terms
		 */
		void set_base(soroban_t new_base)
		{
			score = new_base;
		}

		// using the information already stored in the scoring term variables
		// and the weights supplied to the constructor, calculate the linear combination
		// value and return it
		double calculate_weighted_score()
		{
			return  weights[vdw_weight] * score.vdw +
				weights[vdw_attr_weight] * score.vdw_attr +
				weights[vdw_rep_weight] * score.vdw_rep +
				weights[elec_weight] * score.elec +
				weights[elec_sr_attr_weight] * score.elec_sr_attr +
				weights[elec_lr_attr_weight] * score.elec_lr_attr +
				weights[elec_sr_rep_weight] * score.elec_sr_rep +
				weights[elec_lr_rep_weight] * score.elec_lr_rep +
				weights[hbp_ss_weight] * score.hbp_ss +
				weights[solv_weight] * score.solv +
				weights[sasa_weight] * score.sasa +
				weights[acp_weight] * score.acp;
		}


		// Update the values corresponding to the van der waals components of the score
		void add_vdw(double vdw, double vdw_attr, double vdw_rep)
		{
			this->score.vdw += vdw;
			this->score.vdw_attr += vdw_attr;
			this->score.vdw_rep += vdw_rep;
		}

		// Update the values corresponding to the electrostatics components of the score
		void add_electrostatics(double elec, double elec_sr_attr, double elec_sr_rep, double elec_lr_attr, double elec_lr_rep)
		{
			this->score.elec += elec;
			this->score.elec_sr_attr += elec_sr_attr;
			this->score.elec_sr_rep += elec_sr_rep;
			this->score.elec_lr_attr += elec_lr_attr;
			this->score.elec_lr_rep += elec_lr_rep;
		}

		// Update the values corresponding to the hydrogen bonding potential + disulphide
		void add_hbp(double hbp_ss)
		{
			this->score.hbp_ss += hbp_ss;
		}

		// Update the value corresponding to the acp based desolvation
		void add_acp_based_desolvation(double acp)
		{
			this->score.acp += acp;
		}

		// Update the values corresponding to the solvation components of the score
		void add_solvation(double sasa)
		{
			this->score.sasa += sasa;
		}
		// Update the values corresponding to the lz solvation components of the score
		void add_lz_solvation(double solv)
		{
			this->score.solv += solv;
		}

		/*
		 * Sends a string representation of this object to the output stream
		 */    
		void print(ostream& output_stream)
    {
      print(output_stream, ",");
    }

		void print(ostream& output_stream, string separator)
		{
			output_stream << score.vdw << separator
				<< score.vdw_attr << separator
				<< score.vdw_rep << separator
				<< score.elec << separator
				<< score.elec_sr_attr << separator
				<< score.elec_lr_attr << separator
				<< score.elec_sr_rep << separator
				<< score.elec_lr_rep << separator
				<< score.hbp_ss << separator
				<< score.solv << separator
				<< score.sasa << separator
				<< score.acp;
		}
};

#endif
