#ifndef _ATOM_H_
#define _ATOM_H_

#include <string>
#include <cctype>
#include "charmm.h"

using namespace std;

class atom
{
	public:
	int	anum;   // id of the atom
	string	atype;  // type of the atom
	int	rnum;   // residue id
	double	axyz[3]; // xyz coordinates
	double occupancy; // PDB occupancy field
	double temp_factor; // PDB tempFactor
	string      residue; // name of residue
	string      element_symbol;
	string      chain;// chain id
	double      arad;
        int         INTFA; // interface CA atom
	// These were added as part of SOROBAN, check functionality
	int         INTF;
	double      sas_area;
	int         acp_type; // atom type for Atom contact potential
	int         don_acc; // donor or acceptor
	charmm_param chpar; // charmm values
	double      atomcharge;
	double      atom_solv;//atom solvation parameter
	// true if, during the clash finding process, it is determined that this atom is clashing
	bool isclashing;

	atom(string m_type, string m_residue, string m_chain, int m_rid, double* m_xyz, double m_charge,
         double m_rad, int m_acp_type, charmm_param& m_chpar, string m_element_symbol)
	{
		rnum = m_rid;
		atype = m_type;
		for(int i=0;i<3;i++)
		{
			axyz[i] = m_xyz[i];
		}
		residue = m_residue;
		element_symbol = m_element_symbol;
		arad = m_rad;
		atomcharge = m_charge;
		chpar = m_chpar;
		chain = m_chain;
		INTF = 0;
		acp_type = m_acp_type;
		INTFA=-1;
		isclashing = false;
	}

	atom(int m_num, string m_type, int m_rid, double* m_xyz, string m_residue, string m_chain, double m_solv,
		double m_rad, double m_charge, int m_acp_type, charmm_param& m_p, double m_occupancy, double m_temp_factor,
		string m_element_symbol)
	{
		anum = m_num;
		rnum = m_rid;
		atype = m_type;
		for(int i=0;i<3;i++)
		{
			axyz[i] = m_xyz[i];
		}
		residue = m_residue;
		element_symbol = m_element_symbol;
		chain = m_chain;
		arad = m_rad;
		INTFA=-1;
		INTF = 0;
		// Initialize params added because of SOROBAN
		atom_solv = m_solv;
		acp_type = m_acp_type;
		atomcharge = m_charge;
		chpar = m_p;
		occupancy = m_occupancy;
		temp_factor = m_temp_factor;
		isclashing = false;
	}

	/*
	 * Since we add hydrogens to the PDB's for scoring purposes,
	 * it is useful to have a way of telling if an atom is a hydrogen,
	 * specially when calculating clashes
	 */
	bool is_hydrogen() {
		/* first create a string without the numbers */
		string withoutnumbers;
		for(size_t i = 0; i < atype.length(); i++)
		{
			// add the character if it's not a number
			if(!isdigit(atype[i]))
			{
				withoutnumbers += atype[i];
			}
		}
		// compare it against all known of light atoms coming from hbplus
		return withoutnumbers.compare("H") == 0 ||
			withoutnumbers.compare("HH") == 0 ||
			withoutnumbers.compare("HE") == 0 ||
			withoutnumbers.compare("HD") == 0 ||
			withoutnumbers.compare("HG") == 0 ||
			withoutnumbers.compare("HZ") == 0;
	}

	string tostring() {
		char buf[200];
		sprintf(buf, "ATOM  %5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          ",
		                        anum, atype.c_str(), residue.c_str(), chain.c_str(), rnum, axyz[0], axyz[1], axyz[2],
					                        occupancy, temp_factor);
		if(element_symbol.length() == 1) {
			return string(buf) + " " + element_symbol + "\n";
		} else {
			return string(buf) + element_symbol + "\n";
		}
	}

	atom(){}
	~atom(){}
};

#endif

