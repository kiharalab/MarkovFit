#include "mrf_features.h"

using std::endl;

mrf_one_body_features::mrf_one_body_features( double alpha_z, double beta_y, double gamma_x, double translate_x,
                                             double translate_y, double translate_z,double rmsd, double correlation,
                                             double overlap, string label) : alpha_z(alpha_z), beta_y(beta_y), gamma_x(gamma_x),
translate_x(translate_x), translate_y(translate_y),
translate_z(translate_z), rmsd(rmsd), correlation(correlation),
overlap(overlap), label(label) {
}

void mrf_one_body_features::write(ofstream& data_stream) {
    data_stream << rmsd << "\t" << correlation << "\t" << overlap << "\t"
    << alpha_z << "\t" << beta_y << "\t" << gamma_x << "\t"
    << translate_x << "\t"<< translate_y << "\t"<< translate_z << "\t"
    << label << "\t"
    << endl;
}

void mrf_one_body_features::write_header(ofstream& data_stream) {
    data_stream << "RMSD\tCC\tOverlap\tAlpha_z\tBeta_y\tGamma_x\ttx\tty\ttz\tUnit"
    << endl;
}

mrf_two_body_features::mrf_two_body_features(double alpha_z_left,
                                             double beta_y_left, double gamma_x_left, double translate_x_left,
                                             double translate_y_left, double translate_z_left, double alpha_z_right,
                                             double beta_y_right, double gamma_x_right, double translate_x_right,
                                             double translate_y_right, double translate_z_right, double rmsd,
                                             soroban_score physics_score, int no_clashes, string label_left, string label_right) :
alpha_z_left(alpha_z_left), beta_y_left(beta_y_left),
gamma_x_left(gamma_x_left), translate_x_left(translate_x_left),
translate_y_left(translate_y_left), translate_z_left(translate_z_left),
alpha_z_right(alpha_z_right), beta_y_right(beta_y_right),
gamma_x_right(gamma_x_right), translate_x_right(translate_x_right),
translate_y_right(translate_y_right),
translate_z_right(translate_z_right), rmsd(rmsd),
physics_score(physics_score), no_clashes(no_clashes),
label_left(label_left), label_right(label_right) {
}

void mrf_two_body_features::write(ofstream& data_stream) {
    data_stream << rmsd << "\t";
    data_stream << "\t" << physics_score.calculate_weighted_score() << "\t";
    physics_score.print(data_stream, "\t");
    data_stream << "\t" << no_clashes << "\t" << alpha_z_left << "\t" << beta_y_left << "\t" << gamma_x_left << "\t" << translate_x_left
    << "\t" << translate_y_left << "\t" << translate_z_left << "\t" << alpha_z_right << "\t" << beta_y_right << "\t"
    << gamma_x_right << "\t" << translate_x_right << "\t" << translate_y_right << "\t" << translate_z_right << "\t"
    << label_left << "\t" << label_right << endl;
}

void mrf_two_body_features::write_header(ofstream& data_stream) {
    data_stream << "RMSD\tPhysicsScore\tvdw\tvdw_a\tvdw_r\telec\telec_sr_a\telec_lr_a\telec_sr_r\telec_lr_r\t"
    << "hb\tsolv\tsasa\tacp\tno_clashes\t"
    << "Alpha_z_l\tBeta_y_l\tGamma_x_l\ttx_l\tty_l\ttz_l\tAlpha_z_r\tBeta_y_r\tGamma_x_r\ttx_r\tty_r\ttz_r\t"
    << "Unit_l\tUnit_r" << endl;
}

int mrf_two_body_features::get_no_clashes(){
    return no_clashes;
}

double mrf_two_body_features::get_physics_score(){
    return physics_score.calculate_weighted_score();
}

