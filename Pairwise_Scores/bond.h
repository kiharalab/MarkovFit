#ifndef _BOND_H_
#define _BOND_H_

class bond
{
	public:
	int bstart;
	int bend;
	bond(int mstart, int mend)
	{
		bstart = mstart;
		bend = mend;
	}
	bond(){}
	~bond(){}
};

#endif

