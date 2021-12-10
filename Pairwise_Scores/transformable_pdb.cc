#include "transformable_pdb.h"

void transformable_pdb::scale(double x_scale, double y_scale, double z_scale) {
  for (size_t atom_index = 0; atom_index < this->atoms.size(); atom_index++) {
    this->atoms[atom_index].axyz[0] *= x_scale;
    this->atoms[atom_index].axyz[1] *= y_scale;
    this->atoms[atom_index].axyz[2] *= z_scale;
  }
}

void transformable_pdb::translate(double x_translation, double y_translation, double z_translation) {
  for (size_t atom_index = 0; atom_index < this->atoms.size(); atom_index++) {
    this->atoms[atom_index].axyz[0] += x_translation;
    this->atoms[atom_index].axyz[1] += y_translation;
    this->atoms[atom_index].axyz[2] += z_translation;
  }
}

void transformable_pdb::rotate(double alpha_z, double beta_y, double gamma_x) {
  point_transformation t(alpha_z, beta_y, gamma_x, 0, 0, 0);
  apply_point_transformation(&t);
}

void transformable_pdb::apply_point_transformation( point_transformation* transformation) {
  apply_any_point_transformation(transformation, false);
}

void transformable_pdb::apply_inverse_point_transformation( point_transformation* transformation) {
  apply_any_point_transformation(transformation, true);
}

void transformable_pdb::apply_any_point_transformation(point_transformation* transformation, bool is_inverse) {
  double tmp_x = 0, tmp_y = 0, tmp_z = 0;
  for (size_t atom_index = 0; atom_index < this->atoms.size(); atom_index++) {
    if(is_inverse) {
      transformation->invert_transform(this->atoms[atom_index].axyz[0],
                                       this->atoms[atom_index].axyz[1],
                                       this->atoms[atom_index].axyz[2],
                                       &tmp_x, &tmp_y, &tmp_z);
    } else {
      transformation->transform(this->atoms[atom_index].axyz[0],
                                this->atoms[atom_index].axyz[1],
                                this->atoms[atom_index].axyz[2],
                                &tmp_x, &tmp_y, &tmp_z);
    }
    this->atoms[atom_index].axyz[0] = tmp_x;
    this->atoms[atom_index].axyz[1] = tmp_y;
    this->atoms[atom_index].axyz[2] = tmp_z;
  }
}


void transformable_pdb::apply_point_transformation_sequence(point_transformation_sequence* sequence) {
  double tmp_x = 0, tmp_y = 0, tmp_z = 0;
  for (size_t atom_index = 0; atom_index < this->atoms.size(); atom_index++) {
    sequence->transform(this->atoms[atom_index].axyz[0],
                        this->atoms[atom_index].axyz[1],
                        this->atoms[atom_index].axyz[2],
                        &tmp_x, &tmp_y, &tmp_z);
    this->atoms[atom_index].axyz[0] = tmp_x;
    this->atoms[atom_index].axyz[1] = tmp_y;
    this->atoms[atom_index].axyz[2] = tmp_z;
  }
}

void transformable_pdb::get_center(double & center_x, double & center_y, double & center_z) {
    center_x = 0,center_y = 0, center_z = 0;
    for (size_t atom_index = 0; atom_index < this->atoms.size(); atom_index++) {
        center_x += this->atoms[atom_index].axyz[0];
        center_y += this->atoms[atom_index].axyz[1];
        center_z += this->atoms[atom_index].axyz[2];
    }
    center_x = center_x / this->atoms.size();
    center_y = center_y / this->atoms.size();
    center_z = center_z / this->atoms.size();
}
