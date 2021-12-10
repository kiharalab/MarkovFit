/**
 * PDB header file
 * stores atom coordinates,types and bond information
 * 
 * @author vishwesh venkatraman
 * @date   19/08/2007
 */

#ifndef _PDB_H_
#define _PDB_H_

#include "atom.h"
#include "bond.h"
#include "utils.h"
#include "charmm.h"
#include <vector>
#include <sstream>
#include <map>

// This is according to CAPRI papers (3 angstroms)
#define CLASH_DISTANCE 3.0
// Extracted from the charmm information
#define MAX_ATOM_RADIUS 2.265

class pdb
{
  private:
    // Used to store the centroid obtained by calling calculate_centroid.
    // This is to avoid re-computing it in programs where the centroid is required several times.
    // Every time precompute_centroid is called it will reset these values.
    double centroid_x, centroid_y, centroid_z;
    // Flag read by get_centroid to determine if precompute_centroid has been called before.
    // If it hasn't centroid_? have not been assigned so far Upon construction it's set to false.
    bool is_centroid_precomputed;
	public:
		string pname;
		vector<atom> atoms;
		vector<bond> bonds;

    pdb()
    {
      is_centroid_precomputed = false;
    }

		pdb(string m_name, vector<atom>& m_atoms, vector<bond>& m_bonds)
		{
			pname = m_name;
			atoms = m_atoms;
			bonds = m_bonds;
      is_centroid_precomputed = false;
		}
		
		pdb(string m_name, vector<atom>& m_atoms)
		{
			pname = m_name;
			atoms = m_atoms;
      is_centroid_precomputed = false;
		}

  /*
   * Allocates an empty PDB and merges two PDB's respecting the residue number ordering in each of the PDB's.
   * Originally this process was done by merely concatenating the 2 atom vectors. However, the alignment functions assume
   * that residue numbers are consecutive in each pdb instance (for efficiency purposes).
   * This function makes sure that the order is maintained.
   * The caller must delete the instance returned when finished.
   */
//  static pdb* merge(pdb& first_pdb, pdb& second_pdb);
//  static pdb* reverse_merge(pdb& first_pdb, pdb& second_pdb);

    void merge(pdb& first_pdb, pdb& second_pdb, pdb& both_pdb);
    void reverse_merge(pdb& first_pdb, pdb& second_pdb, pdb& both_pdb);


    // Returns the lowest and highest x, y and z values (in PDB coordinate space). In practice, this is a bounding box enclosing the PDB.

    void get_boundaries(float* min_x, float* min_y, float* min_z, float* max_x, float* max_y, float* max_z);

    /*
     * Compare this.atoms against matched_pdb.atoms and find the c-alpha atoms that exactly match between them.
     * A match is defined as two c-alphas that have the same residue identifier, same 3-letter code and (optionally)
     * the same chain ID.
     *
     * When calling somepdb.get_calphas(other_pdb,...) it is expected that all atoms in somepdb have a correspondence
     * in other_pdb. e.g. if there is a smaller PDB out of the two then the call should be small.get_matching_calphas(bigger,...).
     *
     * Returns as output values copies of the atoms from this instance and also copies from the matched pdb.
     *
     * Note: ignore_chain_id==true => only residue ID and 3-letter code need to be the same for a match to happen.
     */
    void get_matching_calphas(pdb* matched_pdb, bool ignore_chain_id, vector<atom>* atoms_this, vector<atom>* atoms_matched);

    vector<vector<atom> > get_all_matching_calphas(pdb* matched_pdb, vector<atom>* atoms_this);
    vector<vector<atom> > get_all_matching_calphas_multiple(vector<pdb>& aligned_calphas, vector<atom>* atoms_this);
    vector<vector<atom> > get_all_matching_docked_calphas(vector<pdb>& aligned_calphas, vector<atom>* sub_units);

  /*
   * Go through all non-hydrogen atoms and calculate the geometric center (centroid) by averaging the coordinates
   * from each axis, returning the values as output parameters.
   */
  void calculate_centroid(double* x, double* y, double* z);
  // Initializes the centroid_? variables by calling calculate_centroid
  void precompute_centroid();
  // Return the precomputed centroid_? values. If they have not been assigned before, it calls precompute_centroid to do so
  void get_centroid(double* x, double* y, double* z);
};

/**
 * Read PQR informtion from file
 * Parse only "ATOM" information.
 *
 * @see    process_files
 * @params infile name of pqr file
 * @params pdbdata stores pqr coordinates and atom types
 * @params onlycalpha if set to true, the method will only return a collection of C-alpha atoms
 * ignoring all others (by default it ignores this restriction)
 * @return void
 */
void read_protein(string infile, pdb& pdbdata, bool onlycalpha=false);

void write_complex(string pdb_file_name, vector<atom>& RESULT);
void write_complex(string pdb_file_name, vector<vector<atom> >& RESULT);
void write_complex(string fname, int k, vector<atom>& RESULT);
void write_complex(string fname, int k, vector<vector<atom> >& RESULT);

void get_charmm_type(string res, string pdbtype, vector<charmm_type>& CTYPEDATA,
    string& chtype, double& arad, double& charge, int& acptype);

// 0 - none
// 1 - done
// 2 - acceptor
// 3 - both donor and acceptor
/*
Table 1. Hydrogen bond donors and acceptors
Donors
(1) N (main-chain N-H)
(2) Asn OD1, His NE2, His ND1, His CD2(a), His CE1(a), Lys NZ, Asn ND2, Gln NE2 Arg NE, Arg NH1, Arg NH2, Ser OG, Thr OG1, Tyr OH, Trp NE1
Acceptors
(1) O (main-chain C = O)
(2) Asp OD1, Asp OD2, Glu OE1, Glu OE2, His ND1, His CD2a, His CE1(a), Asn OD1, Gln OE1, Gln NE2, Asn ND2, Ser OG, Thr OG1, Tyr OH
(a)We include also the carbon atoms CE1 and CD2 in the prediction of Euler angles from hydrogen bonds to take an ambigous orientation
   of the His imidazole into account.
Meyer et al., J Mol Biol, 1996, 264, 199-210
*/
void set_donor_acceptor(string res, string chtype, int& ac_do);

// based on the atomic desolvation parameters of Abagyan    
double get_adp_parameter(string res, string chtype);

// read charmm atom types from the file "charmm_atom_types.in"
// read charmm parameters from the file "charmm_params.in"
void read_charmm_parameters(vector<charmm_type>& CTYPEDATA, map<string, charmm_param>& CPARAMDATA);

int get_clashing_atoms(vector< vector<atom> >& chains);

int get_no_clashing_atoms(const vector< vector<atom> >& chains);

#endif

