#include "contact.h"

// while looking for the c-alpha it will also set to 1 the INTF flag
// for all the residues that are on the interface.
// INTFA will only identify those residues on the interface and that are C-alpha at the same time
void set_interface_calpha(vector<atom>& X, string chain, string res, int rnum)
{
	for(size_t i=0;i<X.size();i++)
	{
		// all the atoms belonging to the residue defined by chain/res/rnum
		// will be marked as INTF=1
		if(res == X[i].residue && rnum == X[i].rnum && chain == X[i].chain)
		{
			X[i].INTF = 1;
			// and also, to identify an atom as a c-alpha interface one,
			// the INTFA flag is set to 1
			if(X[i].atype == "CA")
			{
				X[i].INTFA = 1;
			}
		}
	}
}

void get_calpha_pos(vector<atom>& X, ANNkd_tree* T)
{
    ANNpoint f = annAllocPt(3);
    double r2 = 100.0;
    for(size_t i=0;i<X.size();i++)
    {
        f[0] = X[i].axyz[0];
        f[1] = X[i].axyz[1];
        f[2] = X[i].axyz[2];
        ANNidxArray nnIdx = NULL;
        ANNdistArray dists = NULL;

        int nnode = T->annkFRSearch(f, r2, 0, nnIdx, dists, 0.); 
        // set the CA atom of the residue to which this atom belongs            
        if(nnode > 0)
	{
            set_interface_calpha(X, X[i].chain, X[i].residue, X[i].rnum);
	}
    }
}

void get_interface_residues(vector<atom>& R, vector<atom>& L)
{
    // get KDTree for R and L
    // get KDTree for receptor
    int N = (int) R.size();
    ANNpointArray RdataPts;
    RdataPts = annAllocPts(N, 3);

    for(int i=0; i<N; i++)
    {
        for(int j=0;j<3;j++)
            RdataPts[i][j] = R[i].axyz[j];
    }
    ANNkd_tree *TR = new ANNkd_tree(RdataPts, N, 3);

    N = (int) L.size();
    ANNpointArray LdataPts;
    LdataPts = annAllocPts(N, 3);

    for(int i=0; i<N; i++)
    {
        for(int j=0;j<3;j++)
	{
            LdataPts[i][j] = L[i].axyz[j];
	}
    }
    ANNkd_tree *TL = new ANNkd_tree(LdataPts, N, 3);

    get_calpha_pos(R, TL);
    get_calpha_pos(L, TR);

    // kdtree deletion
    delete TR;
    annDeallocPts(RdataPts);
    annClose();    
    // kdtree deletion
    delete TL;
    annDeallocPts(LdataPts);
    annClose();
}

