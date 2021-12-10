/*
 * Derivation from the original pdb class that allows the modification
 * of atomic coordinates using scaling, rotation and translation operations
 */
#ifndef _TRANSFORMABLE_PDB_H_
#define _TRANSFORMABLE_PDB_H_

#include "point_transformation.h"
#include "pdb.h"

class transformable_pdb : public pdb
{
 public:
  transformable_pdb() : pdb() {}
  transformable_pdb(string m_name, vector<atom>& m_atoms) :
	    pdb(m_name, m_atoms) {}
  // Scale all atomic coordinates in the parent's atoms collection,
  // replacing them by the transformed coordinates
  void scale(double x_scale, double y_scale, double z_scale);
  // Adds the specified magnitude to every atomic coordinate, overwritting the previous values.
  void translate(double x_translation, double y_translation, double z_translation);
  // Apply a point_transformation without translation, only rotation.
  void rotate(double alpha_z, double beta_y, double gamma_x);
  // Apply transformation.transform to overwrite all atom coordinates.
  // A method to apply the inverse transformation is also provided.
  void apply_point_transformation(point_transformation* transformation);
  void apply_inverse_point_transformation(point_transformation* transformation);
  // Apply a sequence of transformations calling
  // point_transformation_sequence.transform on every atom
  void apply_point_transformation_sequence(point_transformation_sequence* sequence);
  void get_center(double & center_x, double & center_y, double & center_z);

 private:
  // Auxiliary used by both apply_transformation methods. Normal
  // transformation when the second parameter is false, and inverse when true.
  //  (just to keep almost the same code in the same place).
  void apply_any_point_transformation(point_transformation* transformation, bool is_inverse);

};

#endif
