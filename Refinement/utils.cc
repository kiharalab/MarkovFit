#include "utils.h"

#include <sstream>

using std::string;
using std::stringstream;
using std::vector;

string trimmed(string original, int num_chars)
{
	size_t first_char_pos = original.find_first_not_of(' ');
	size_t last_char_pos = original.find_last_not_of(' ');
	if(first_char_pos == string::npos) {
		return "";
	} else {
		size_t characters_kept = last_char_pos - first_char_pos + 1;
		return original.substr(first_char_pos, characters_kept);
	}
}

string int_to_string(int t)
{
    stringstream ss;
    ss << t;
    return ss.str();
}

vector<string> split(string base, char separator)
{
	vector<string> split_string;
	while(!base.empty())
	{
		// get the following element first
		size_t separator_pos = base.find_first_of(separator);
		string new_element = base.substr(0, separator_pos);
		split_string.push_back(new_element);
		// and then remove the chain and the separator before the next iteration starts
		base = separator_pos ==  string::npos ? "" : base.substr(separator_pos + 1);
	}
	return split_string;
}
/*
 * Calculate euclidean distance between points. Each double* points to a 3D array that
 * contains the x,y,z coordinates, respectively.
 */
double get_distance(double *p, double *q)
{
	double d = sqrt(get_squared_distance(p, q));
	return d;
}

/*
 * Calculate euclidean distance between points and then return that distance squared.
 * This is actually the same calculation for the euclidean distance although the square root is not calculated
 * at the end of the function
 * Each double* points to a 3D array that
 * contains the x,y,z coordinates, respectively.
 */
double get_squared_distance(double *p, double *q)
{
	double f0 = p[0] - q[0];
	double f1 = p[1] - q[1];
	double f2 = p[2] - q[2];

	f0 *= f0;f1 *= f1; f2 *= f2;
	double d = f0 + f1 + f2;
	return d;
}
