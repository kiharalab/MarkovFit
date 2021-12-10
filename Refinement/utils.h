#ifndef _UTILS_H_
#define _UTILS_H_

#include <string>
#include <vector>
#include <sstream>
#include <cmath>

using std::string;
using std::vector;


string int_to_string(int t);
string trimmed(string original, int num_chars);
vector<string> split(string base, char separator);
/*
 * Calculate euclidean distance between points. Each double* points to a 3D array that
 * contains the x,y,z coordinates, respectively.
 */
double get_distance(double *p, double *q);
/*
 * Calculate euclidean distance between points and then returned that distance squared.
 * This is actually the same calculation for the euclidean distance although the square root is not calculated
 * at the end of the function
 * Each double* points to a 3D array that
 * contains the x,y,z coordinates, respectively.
 */
double get_squared_distance(double *p, double *q);

#endif
