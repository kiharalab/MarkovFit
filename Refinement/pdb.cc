#include "pdb.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <climits>

#include <ANN/ANN.h>

using namespace std;

// Hydrogen
#define H_RAD 1.20
// carbon
#define C_RAD = 1.70
// nitrogen
#define N_RAD = 1.65
// oxygen
#define O_RAD = 1.60
// phosphorous
#define P_RAD = 1.90
// sulphur
#define S_RAD = 1.80
// default
#define DEF_RAD = 1.70

int kBlankChars;

/*** Methods added to the pdb class ***/

void pdb::merge(pdb& first_pdb, pdb& second_pdb, pdb& merged)
{
    //    pdb* merged = new pdb();
   
    for (size_t i = 0; i < first_pdb.atoms.size(); i++) {
        atom candidate = first_pdb.atoms[i];
        merged.atoms.push_back(candidate);
    }
    for (size_t i = 0; i < second_pdb.atoms.size(); i++) {
        atom candidate = second_pdb.atoms[i];
        merged.atoms.push_back(candidate);
    }
    
    //    return merged;
}

void pdb::reverse_merge(pdb& first_pdb, pdb& second_pdb, pdb& merged)
{
    for (size_t i = 0; i < second_pdb.atoms.size(); i++) {
        atom candidate = second_pdb.atoms[i];
        merged.atoms.push_back(candidate);
    }
    
    for (size_t i = 0; i < first_pdb.atoms.size(); i++) {
        atom candidate = first_pdb.atoms[i];
        merged.atoms.push_back(candidate);
    }
}

void pdb::get_boundaries(float* min_x, float* min_y, float* min_z, float* max_x, float* max_y, float* max_z)
{
    
    (*min_x) = (*min_y) = (*min_z) = LONG_MAX;
    (*max_x) = (*max_y) = (*max_z) = LONG_MIN;
    
    float x, y, z;
    
    for(vector<atom>::iterator it = atoms.begin(); it != atoms.end(); it++) {
        atom current_atom = *it;
        
        x = current_atom.axyz[0];
        y = current_atom.axyz[1];
        z = current_atom.axyz[2];
        
        if (x < (*min_x)) (*min_x) = x;
        if (y < (*min_y)) (*min_y) = y;
        if (z < (*min_z)) (*min_z) = z;
        
        if (x > (*max_x)) (*max_x) = x;
        if (y > (*max_y)) (*max_y) = y;
        if (z > (*max_z)) (*max_z) = z;
    }
}


void pdb::get_matching_calphas(pdb* matched_pdb, bool ignore_chain_id, vector<atom>* atoms_this, vector<atom>* atoms_matched) {
    // get all c-alphas from this
    for(size_t this_index = 0; this_index < atoms.size(); this_index++) {
        atom current_atom = atoms[this_index];
        if(!current_atom.atype.compare("CA")) {
#ifdef DEBUG
            cerr << "Source CA: " << current_atom.residue << " " << current_atom.rnum << endl;
#endif
            atoms_this->push_back(current_atom);
        }
    }
    // find a match for each and every c-alpha just retrieved
    size_t matched_index = 0;
    size_t calpha_index = 0;
    for(; calpha_index < atoms_this->size(); calpha_index++) {
        atom current_calpha = (*atoms_this)[calpha_index];
        for(; matched_index < matched_pdb->atoms.size(); matched_index++) {
            atom matched_atom = matched_pdb->atoms[matched_index];
#ifdef DEBUG
            cerr << "Comparison: " << current_calpha.residue << " " << current_calpha.rnum << " "
            << matched_atom.residue << " " << matched_atom.rnum << endl;
#endif
            if((!matched_atom.atype.compare("CA"))
               && matched_atom.rnum == current_calpha.rnum
               && (!matched_atom.residue.compare(current_calpha.residue))
               && (ignore_chain_id || (!matched_atom.chain.compare(current_calpha.chain)))) {
                //match was found,allow outer loop to go to next atom and advanced index in this one for next comparison to be ready
#ifdef DEBUG
                cerr << "Match: " << current_calpha.residue << " " << current_calpha.rnum << endl;
#endif
                atoms_matched->push_back(matched_atom);
                matched_index++;
                break;
            }
        }
    }
#ifdef DEBUG
    cerr << "Current c-alpha index " << calpha_index << " Total " << atoms_this->size() << endl;
#endif
    // if any of the atoms was not correctly matched then return empty collections (because partial matches are not expected)
    if(calpha_index < atoms_this->size()) {
        atoms_this->clear();
        atoms_matched->clear();
    }
}

vector<vector<atom> > pdb::get_all_matching_calphas(pdb* matched_pdb, vector<atom>* atoms_this) {
    // get all c-alphas from this
    for(size_t this_index = 0; this_index < atoms.size(); this_index++) {
        atom current_atom = atoms[this_index];
        if(!current_atom.atype.compare("CA")) {
            //cerr << "Source CA: " << current_atom.residue << " " << current_atom.rnum << endl;
            atoms_this->push_back(current_atom);
        }
    }
    // find a match for each and every c-alpha just retrieved
    size_t matched_index = 0;
    size_t calpha_index = 0;
    //cout<< "Before anything, size = "<< atoms_this->size() <<" ,  " << (*atoms_this)[0].chain << endl;
    
    vector<atom> atoms_matched;
    vector<vector<atom> > all_atoms_matched;
    size_t count = 0;
    
    for(; calpha_index < atoms_this->size(); calpha_index++) {
        atom current_calpha = (*atoms_this)[calpha_index];
        for(; matched_index < matched_pdb->atoms.size(); matched_index++) {
            atom matched_atom = matched_pdb->atoms[matched_index];
            
            if((!matched_atom.atype.compare("CA"))
               && matched_atom.rnum == current_calpha.rnum
               && (!matched_atom.residue.compare(current_calpha.residue))) {
                count=0;
                for (size_t i=0; i<atoms_this->size(); i++) {
                    atom sub_atom = (*atoms_this)[calpha_index+i];
                    atom struct_atom = matched_pdb->atoms[matched_index+i];
                    if((!struct_atom.atype.compare("CA"))
                       && sub_atom.rnum == struct_atom.rnum
                       && (!sub_atom.residue.compare(struct_atom.residue))) {
                        atoms_matched.push_back(struct_atom);
                        count++;
                        if(atoms_matched.size() == atoms_this->size()) {
                            all_atoms_matched.push_back(atoms_matched);
                            cout<< sub_atom.chain <<" is matched with " << struct_atom.chain << endl;
                            cout<< "matched index = " << matched_index <<" , total = " << matched_pdb->atoms.size() << endl;
                            return all_atoms_matched;
                        }
                        
                    } else {
                        atoms_matched.clear();
                        calpha_index=-1;
                        matched_index++;
                        break;
                    }
                    
                }
                break;
            }
        }
    }
    
    return all_atoms_matched;
}

vector<vector<atom> > pdb::get_all_matching_calphas_multiple(vector<pdb>& aligned_calphas, vector<atom>* sub_units) {
    // get all c-alphas from subunits
    for(size_t sub_index = 0; sub_index < atoms.size(); sub_index++) {
        atom current_atom = atoms[sub_index];
        if(!current_atom.atype.compare("CA")) {
            sub_units->push_back(current_atom);
        }
    }
    
    vector<atom> matched_atoms;
    vector<vector<atom> > all_matchings;
    
    for(size_t i=0; i< aligned_calphas.size(); i++) {
        pdb current_pdb = aligned_calphas[i];
        if(current_pdb.atoms.size() != sub_units->size()) {
            continue;
        }
        
        for(size_t k=0; k<current_pdb.atoms.size(); k++) {
            atom sub_atom = (*sub_units)[k];
            atom all_atom = current_pdb.atoms[k];
            
            if(
               (!all_atom.atype.compare("CA"))
               //&& sub_atom.rnum == all_atom.rnum
               &&
               (!sub_atom.residue.compare(all_atom.residue))) {
                matched_atoms.push_back(all_atom);
            }
        }
        
        if(matched_atoms.size() == current_pdb.atoms.size()) {
            all_matchings.push_back(matched_atoms);
        }
        matched_atoms.clear();
    }
    return all_matchings;
}

vector<vector<atom> > pdb::get_all_matching_docked_calphas(vector<pdb>& aligned_calphas, vector<atom>* sub_units) {
    // get all c-alphas from subunits
    sub_units->clear();
    for(size_t sub_index = 0; sub_index < atoms.size(); sub_index++) {
        atom current_atom = atoms[sub_index];
        if(!current_atom.atype.compare("CA")) {
            sub_units->push_back(current_atom);
        }
    }
    
    vector<atom> matched_atoms;
    vector<atom> reversed_matched_atoms;
    vector<vector<atom> > all_matchings;
    
    for(size_t i=0; i< aligned_calphas.size()-1; i++) {
        for(size_t j=i+1; j< aligned_calphas.size(); j++) {
            pdb both_pdb, reversed_both_pdb;
            both_pdb.merge(aligned_calphas[i], aligned_calphas[j],both_pdb);
            reversed_both_pdb.reverse_merge(aligned_calphas[i], aligned_calphas[j],reversed_both_pdb);

            if(both_pdb.atoms.size() != sub_units->size()) {
                continue;
            }
            
            for(size_t k=0; k<both_pdb.atoms.size(); k++) {
                atom sub_atom = (*sub_units)[k];
                atom all_atom = both_pdb.atoms[k];
                atom r_all_atom = reversed_both_pdb.atoms[k];
                if(
                   (!all_atom.atype.compare("CA"))
                   &&
                   (!sub_atom.residue.compare(all_atom.residue))) {

                    matched_atoms.push_back(all_atom);
                }
                
                if(
                   (!r_all_atom.atype.compare("CA"))
                   &&
                   (!sub_atom.residue.compare(r_all_atom.residue))) {
                    reversed_matched_atoms.push_back(r_all_atom);
                }
            }

            if(matched_atoms.size() == both_pdb.atoms.size()) {
                all_matchings.push_back(matched_atoms);
            }
            
            if(reversed_matched_atoms.size() == reversed_both_pdb.atoms.size()) {
                all_matchings.push_back(reversed_matched_atoms);
            }
            matched_atoms.clear();
            reversed_matched_atoms.clear();
        }
    }

    vector<atom>().swap( matched_atoms );
    vector<atom>().swap( reversed_matched_atoms );
    return all_matchings;
    
}

void pdb::calculate_centroid(double* x, double* y, double* z) {
    double partial_sum_x = 0.0f;
    double partial_sum_y = 0.0f;
    double partial_sum_z = 0.0f;
    
    for(size_t atom_index = 0; atom_index < this->atoms.size(); atom_index++) {
        atom current_atom = this->atoms[atom_index];
        partial_sum_x += current_atom.axyz[0];
        partial_sum_y += current_atom.axyz[1];
        partial_sum_z += current_atom.axyz[2];
    }
    
    *x = partial_sum_x / static_cast<double>(this->atoms.size());
    *y = partial_sum_y / static_cast<double>(this->atoms.size());
    *z = partial_sum_z / static_cast<double>(this->atoms.size());
}

void pdb::precompute_centroid() {
    calculate_centroid(&centroid_x, &centroid_y, &centroid_z);
}

void pdb::get_centroid(double* x, double* y, double* z) {
    if (!is_centroid_precomputed) {
        precompute_centroid();
    }
    *x = centroid_x;
    *y = centroid_y;
    *z = centroid_z;
}

/*** Functions not contained by any particular class ***/
void print_pdb(FILE* pdb_file, vector<atom>& L)
{
    // This used to add TER when we switch to a new chain ID
    string last_chain = "";
    for(size_t i=0;i<L.size();i++)
    {
        if(!last_chain.empty() && last_chain.compare(L[i].chain)) // this returns non-zero if they are different, in which case we add TER
        {
            fprintf(pdb_file, "TER\n");
        }
        
        fprintf(pdb_file, "ATOM  %5d  %-4s%3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          ",
                L[i].anum, L[i].atype.c_str(), L[i].residue.c_str(), L[i].chain.c_str(), L[i].rnum, L[i].axyz[0], L[i].axyz[1], L[i].axyz[2],
                L[i].occupancy, L[i].temp_factor);
        // print an empty space first
        if(L[i].element_symbol.length() == 1) {
            fprintf(pdb_file, " %s\n", L[i].element_symbol.c_str());
        } else {
            fprintf(pdb_file, "%s\n", L[i].element_symbol.c_str());
        }
        last_chain = L[i].chain;
    }
}

void write_complex(string fname, int k, vector<atom>& RESULT)
{
    vector<vector<atom> > wrapper;
    wrapper.push_back(RESULT);
    write_complex(fname, k, wrapper);
}

void write_complex(string fname, int k, vector<vector<atom> >& RESULT)
{
    char numstring[12];
    sprintf(numstring, "%03d", k + 1); /* +1 to make it more readable (not zero-indexed) */
    string pdb_file_name = fname + "_" + string(numstring) + ".pdb";
    write_complex(pdb_file_name, RESULT);
}

void write_complex(string pdb_file_name, vector<atom>& RESULT)
{
    vector<vector<atom> > wrapper;
    wrapper.push_back(RESULT);
    write_complex(pdb_file_name, wrapper);
}

void write_complex(string pdb_file_name, vector<vector<atom> >& RESULT)
{
    // write pdb
    FILE *pdb_file;
    if(pdb_file_name.length() == 0)
    {
        /* if it wasn't provided sent it to stdout */
        pdb_file = stdout;
    }
    else if((pdb_file = fopen( pdb_file_name.c_str(), "w")) == NULL)
    {
        printf( "This file could not be opened.\nDying\n\n" ) ;
        exit(EXIT_FAILURE);
    }
    for(size_t i=0;i<RESULT.size();i++)
    {
        vector<atom> F = RESULT[i];
        print_pdb(pdb_file, F);
        if(i<(RESULT.size()-1))
            fputs("TER\n", pdb_file);
    }
    fclose(pdb_file);
}


void get_charmm_type(string res, string pdbtype, vector<charmm_type>& CTYPEDATA,
                     string& chtype, double& arad, double& charge, int& acptype)
{
    bool found = false;
    for(size_t i=0;i<CTYPEDATA.size();i++)
    {
        if(CTYPEDATA[i].RES == res && CTYPEDATA[i].PDBTYPE == pdbtype)
        {
            chtype = CTYPEDATA[i].CHARMMTYPE;
            arad = CTYPEDATA[i].atomrad;
            charge = CTYPEDATA[i].charge;
            acptype = CTYPEDATA[i].acptype;
            found = true;
            break;
        }
    }
    if(!found)
    {
        cerr << "charmm types not available: " << res << ", " << pdbtype << endl;
        chtype = "NULL";
        arad = 1.9;
        charge = 0.0;
        acptype = 0;
    }
}


// 0 - none
// 1 - done
// 2 - acceptor
// 3 - both donor and acceptor
/*
 Table 1. Hydrogen bond donors and acceptors
 Donors
 (1) N (main-chain N-H)
 (2) Asn OD1, His NE2, His ND1, His CD2(a), His CE1(a), Lys NZ, Asn ND2, Gln NE2 Arg NE, Arg NH1, Arg NH2, Ser OG, Thr OG1,
 Tyr OH, Trp NE1
 Acceptors
 (1) O (main-chain C = O)
 (2) Asp OD1, Asp OD2, Glu OE1, Glu OE2, His ND1, His CD2a, His CE1(a), Asn OD1, Gln OE1, Gln NE2, Asn ND2, Ser OG, Thr
 OG1, Tyr OH
 (a)We include also the carbon atoms CE1 and CD2 in the prediction of Euler angles from hydrogen bonds to take an ambigous
 orientation of the His imidazole into account.
 
 Meyer et al., J Mol Biol, 1996, 264, 199-210
 */
void set_donor_acceptor(string res, string chtype, int& ac_do)
{
    if(res == "HIS" && (chtype == "ND1" || chtype == "CD2" || chtype == "CE1"))
    {
        ac_do = 3;
        return;
    }
    if(res == "HIS" && (chtype == "NE2"))
    {
        ac_do = 1;
        return;
    }
    if(res == "THR" && chtype == "OG1")
    {
        ac_do = 3;
        return;
    }
    if(res == "TYR" && chtype == "OH")
    {
        ac_do = 3;
        return;
    }
    if(res == "TRP" && chtype == "NE1")
    {
        ac_do = 1;
        return;
    }
    if(res == "ASN" && (chtype == "ND2" || chtype == "OD1"))
    {
        ac_do = 3;
        return;
    }
    if(res == "LYS" && chtype == "NZ")
    {
        ac_do = 1;
        return;
    }
    if(res == "SER" && chtype == "OG")
    {
        ac_do = 3;
        return;
    }
    if(res == "GLN" && chtype == "NE2")
    {
        ac_do = 3;
        return;
    }
    if(res == "GLN" && chtype == "OE1")
    {
        ac_do = 2;
        return;
    }
    if(res == "ARG" && (chtype == "NH1" || chtype == "NH2" || chtype == "NE"))
    {
        ac_do = 1;
        return;
    }
    if(chtype == "N") // backbone amide
    {
        ac_do = 1;
        return;
    }
    
    // acceptor
    if(res == "ASP" && (chtype == "OD1" || chtype == "OD2"))
    {
        ac_do = 2;
        return;
    }
    if(res == "GLU" && (chtype == "OE1" || chtype == "OE2"))
    {
        ac_do = 2;
        return;
    }
    if(chtype == "O")
    {
        ac_do = 2;
        return;
    }
    ac_do = 0;
    return;
}

// based on the atomic desolvation parameters of Abagyan    
double get_adp_parameter(string res, string chtype)
{
    // aromatic carbon
    if(chtype == "CR" || chtype == "CR1E")
    {
        return 110.80;
    }
    // aliphatic carbon
    if(chtype == "CT")
    {
        return 19.18;
    }
    // hydroxyl oxygen
    if(chtype == "OH1")
    {
        return -42.55;
    }
    // O- in Glu, ASP
    if(chtype.at(0) == 'O')
    {
        if((res.compare("ASP") == 0) || (res.compare("GLU") == 0))
            return -68.77;
    }
    // carbonyl oxygen
    if(chtype == "O")
    {
        return -31.28;
    }
    // NE1, NE2 in Arg+
    if(chtype == "NE1" || chtype == "NE2")
    {
        if(res.compare("ARG") == 0)
        {
            return -62.56;
        }
    }
    // NZ in Lys+
    if(chtype == "NH3")
    {
        if(res == "LYS")
        {
            return -126.04;
        }
    }
    // uncharged N
    if(chtype == "N")
    {
        return -39.10;
    }
    // S in SH
    if(chtype == "SH1E")
    {
        return 25.76;
    }
    // S in Met or S-S
    if(chtype == "S")
    {
        return 5.06;
    }
    
    //cerr << "Warning: Atomic solvation parameter not found for " << chtype << ":" << res << endl;
    
    return 0.;
}

void read_protein(string infile, pdb& pdbdata, bool onlycalpha)
{
    vector<charmm_type> CTYPEDATA;
    map<string, charmm_param> CPARAMDATA;
    read_charmm_parameters(CTYPEDATA, CPARAMDATA);
    
    ifstream ifs;
    ifs.open(infile.c_str());
    if(!ifs)
    {
        cerr << "Could not open file: " << infile << endl;
        exit(EXIT_FAILURE);
    }
    
    string line;
    
    int atom_id, res_id;
    string res_name, atom_type, chain_id, element_symbol;
    double d_xyz[3];
    
    double solv_param = 0.;
    double occupancy = 0;
    double temp_factor = 0;
    
    vector<atom> vct_atom;
    
    while(getline(ifs, line))
    {
        if((line.substr(0, 4)).compare("ATOM") == 0)
        {
            atom_id = atoi((line.substr(6, 5)).c_str());
            atom_type = trimmed(line.substr(12, 4), kBlankChars); // hydrogen atoms have 4 letter identification
            res_name = trimmed(line.substr(17, 3), kBlankChars);
            element_symbol = line.length() < 78 ? "" : trimmed(line.substr(76, 2), kBlankChars);
            res_id = atoi((line.substr(22, 4)).c_str());
            chain_id = trimmed(line.substr(21, 1), kBlankChars);
            
            
            if(onlycalpha && atom_type.compare("CA") != 0) // then the current atom is not a c-alpha (we ignore it)
            {
                continue;
            }
            
            d_xyz[0] = atof((line.substr(30, 8)).c_str());
            d_xyz[1] = atof((line.substr(38, 8)).c_str());
            d_xyz[2] = atof((line.substr(46, 8)).c_str());
            
            if(line.length() > 54) // i.e. we don't have occupancy or tempFactor info
            {
                occupancy = atof((line.substr(54, 6)).c_str());
                temp_factor = atof((line.substr(60, 6)).c_str());
            }
            
            // hack for HSE residue type
            if(res_name.compare("HSE") == 0) res_name = "HIS";
            // hack for MSE residue type
            if(res_name.compare("MSE") == 0) res_name = "MET";
            // get corresponding charmm atom type
            string chtype = "";
            double arad, achg;
            int acptype;
            get_charmm_type(res_name, atom_type, CTYPEDATA, chtype, arad, achg, acptype);
            
            charmm_param nc;
            if(chtype.compare("NULL") != 0)
            {
                // lookup chtype
                nc = CPARAMDATA[chtype];
                
                if(chtype.at(0) != 'H')
                {
                    //solv_param = get_solvation_parameter(res_name, chtype);
                    solv_param = get_adp_parameter(res_name, chtype);
                }
            }
            else
            {
                solv_param = 0.;
            }
            
            // create atom object
            atom n_atom(atom_id, atom_type, res_id, d_xyz, res_name, chain_id, solv_param, arad, achg, acptype, nc,
                        occupancy, temp_factor, element_symbol);
            
            if(atom_type.at(0) != 'H')
            {
                // assign acceptor/donor, hybridization states
                int ac_do;
                set_donor_acceptor(res_name, atom_type, ac_do);
                n_atom.don_acc = ac_do;
            }
            
            vct_atom.push_back(n_atom);
        }
    }// end while
    
    ifs.close();
    
    string fname = infile;
    
    pdbdata.atoms = vct_atom;
    pdbdata.pname = fname;
}

// read charmm atom types from the file "charmm_atom_types.in"
// read charmm parameters from the file "charmm_params.in"
void read_charmm_parameters(vector<charmm_type>& CTYPEDATA, map<string, charmm_param>& CPARAMDATA)
{
    ifstream ifs;
    ifs.open(CHARMM_ATOM_TYPES.c_str());
    if(!ifs)
    {
        cerr << "Could not open file: " << CHARMM_ATOM_TYPES << endl;
        exit(EXIT_FAILURE);
    }
    
    string line;
    while(getline(ifs, line))
    {
        if(line.find("#") != string::npos)
        {
            continue;
        }
        istringstream iss(line);
        string s1, s2, s3;
        double chg, rad;
        int acptype;
        // residue, pdb atom type, charmm atom type
        iss >> s1 >> s2 >> s3 >> chg >> rad >> acptype;
        if(s3 == "-")
            s3 = "NULL";
        charmm_type nc(s1, s2, s3, chg, rad, acptype);
        CTYPEDATA.push_back(nc);
    }
    
    ifs.close();
    
    // read params
    ifstream ifs2;
    ifs2.open(CHARMM_ATOM_PARAMS.c_str());
    if(!ifs2)
    {
        cerr << "Could not open file: " << CHARMM_ATOM_PARAMS << endl;
        exit(EXIT_FAILURE);
    }
    
    while(getline(ifs2, line))
    {
        if(line.find("#") != string::npos)
        {
            continue;
        }
        istringstream iss(line);
        string str[8];
        // residue, pdb atom type, charmm atom type
        for(int i=0;i<8;i++) {
            iss >> str[i];
            //cout<<i<<" "<< str[i] <<" ";
        }
        //cout<<"\n";
        charmm_param nc(atof(str[1].c_str()), atof(str[2].c_str()), atof(str[3].c_str()),
                        atof(str[4].c_str()), atof(str[5].c_str()), atof(str[6].c_str()));
        CPARAMDATA[str[0]] = nc;
    }
    
    ifs2.close();
}

/*
 * Less efficient version of the clash counting method whose sole purpose is debugging
 * Thus, it is implemented in a simple way (which makes it less efficient).
 * The objective is to obtain a vector of the atoms that are considered to be clashing
 */
int get_clashing_atoms(vector< vector<atom> >& chains)
{
    vector<atom> clashing_atoms;
//    int no_clashes = 0;
    vector<atom> cumulative_atoms;
    for(size_t chain_index = 0; chain_index < chains.size(); chain_index++)
    {
        vector<atom> current_chain = chains[chain_index];
        /* compare the distances of the atoms in the current chain with all the accumulated ones */
        for(size_t current_chain_index = 0; current_chain_index < current_chain.size(); current_chain_index++)
        {
            atom tested_atom = current_chain[current_chain_index];
            /* ignore hydrogens */
            if(tested_atom.is_hydrogen())
                continue;
            for(size_t cumulative_index = 0; cumulative_index < cumulative_atoms.size(); cumulative_index++)
            {
                atom other_atom = cumulative_atoms[cumulative_index];
                if(other_atom.is_hydrogen())
                    continue;
                /* if the distance is less that the threshold we add it to the list and check the next one*/
                if(get_squared_distance(tested_atom.axyz, other_atom.axyz) < CLASH_DISTANCE * CLASH_DISTANCE)
                {
                    if(!current_chain[current_chain_index].isclashing)
                    {
                        clashing_atoms.push_back(current_chain[current_chain_index]);
                    }
                    /* also check if the other atom has been counted before. If it hasn't, do so */
                    if(!other_atom.isclashing)
                    {
                        clashing_atoms.push_back(cumulative_atoms[cumulative_index]);
                    }
                    /* set the clashing flag for both of them */
                    current_chain[current_chain_index].isclashing = true;
                    cumulative_atoms[cumulative_index].isclashing = true;
                }
            }
        }
        // add the current protein to the current list of atoms
        cumulative_atoms.reserve(cumulative_atoms.size() + current_chain.size());
        cumulative_atoms.insert(cumulative_atoms.begin() + cumulative_atoms.size(),
                                current_chain.begin(), current_chain.end());
    }
//    return no_clashes;
    return clashing_atoms.size();

}

/*
 * Less efficient version of the clash counting method whose sole purpose is debugging
 * Thus, it is implemented in a simple way (which makes it less efficient).
 * The objective is to obtain a vector of the atoms that are considered to be clashing
 */
int get_no_clashing_atoms(const vector< vector<atom> >& chains)
{
    float clash_dist = 3;
    int no_clashes = 0;

    //Find atoms in chains 1 &2
    vector<atom> chain1, chain2;
    for (size_t chain1_index = 0; chain1_index < chains[0].size(); chain1_index++) {
        if (chains[0][chain1_index].atype.at(0) == 'H') {
            continue;
        }
        chain1.push_back(chains[0][chain1_index]);
    }

    for (size_t chain2_index = 0; chain2_index < chains[1].size(); chain2_index++) {
        if (chains[1][chain2_index].atype.at(0) == 'H') { //tested_atom.is_hydrogen()
            continue;
        }
        chain2.push_back(chains[1][chain2_index]);
    }

    vector<atom> residue1;
    int res_no;

    /* compare the distances of the atoms in the current chain with all the accumulated ones */
    for(size_t chain1_index = 0; chain1_index < chain1.size(); chain1_index++)
    {
        residue1.clear();
        res_no=chain1[chain1_index].rnum;
        while (chain1[chain1_index].rnum == res_no) {
            residue1.push_back(chain1[chain1_index]);
            chain1_index++;
        }
        chain1_index--;

        for(size_t chain2_index = 0; chain2_index < chain2.size(); chain2_index++)
        {
            atom other_atom = chain2[chain2_index];

            for (size_t res1_index = 0; res1_index < residue1.size(); res1_index++) {
                atom tested_atom = residue1[res1_index];

                /* if the distance is less that the threshold we add it to the list and check the next one*/
                if (get_squared_distance(tested_atom.axyz, other_atom.axyz) < clash_dist * clash_dist) {
                    no_clashes++;
                    res_no=chain2[chain2_index].rnum;
                    while (chain2[chain2_index].rnum == res_no) {
                        chain2_index++;
                    }
                    chain2_index--;
                    break;
                }//end if
            }//end res1
        }//end chain2
    }//end chain1

    vector<atom>().swap( chain1 );
    vector<atom>().swap( chain2 );
    vector<atom>().swap( residue1 );

    return no_clashes;
}
