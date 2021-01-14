#include <panglos.hpp>

// Energy Grid from NJOY2016, to use in default constructor, in units of eV
std::vector<double> default_energy_grid = {
  1.e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5, 1.e-4, 1.26e-4, 1.6e-4,
  2.0e-4, .000253, .000297, .000350, .00042, .000506, .000615, .00075e0,
  .00087, .001012, .00123, .0015, .0018, .00203, .002277, .0026, .003, .0035,
  .004048, .0045, .005, .0056, .006325, .0072, .0081, .009108, .01, .01063,
  .0115, .012397, .0133, .01417, .015, .016192, .0182, .0199, .020493, .0215,
  .0228, .0253, .028, .030613, .0338, .0365, .0395, .042757, .0465, .050,
  .056925, .0625, .069, .075, .081972, .09, .096, .1035, .111573, .120, .128,
  .1355, .145728, .160, .172, .184437, .20, .2277, .2510392, .2705304,
  .2907501, .3011332, .3206421, .3576813, .39, .4170351, .45, .5032575, .56,
  .625, .70, .78, .86, .95, 1.05, 1.16, 1.28, 1.42, 1.55, 1.70, 1.855, 2.02,
  2.18, 2.36, 2.59, 2.855, 3.12, 3.42, 3.75, 4.07, 4.46, 4.90, 5.35, 5.85,
  6.40, 7.00, 7.65, 8.40, 9.15, 9.85, 10.00
};

// CDF grid used by Pavlou and Ji in their first paper [1].
std::vector<double> default_cdf_grid = {
  0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3,
  0.325, 0.35, 0.375, 0.4, 0.425, 0.45, 0.475, 0.5, 0.525, 0.55, 0.575, 0.6,
  0.625, 0.65, 0.675, 0.7, 0.725, 0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9,
  0.925, 0.95, 0.975
};

Panglos::Panglos(): beta_max_(20.),
                    temperature_(),
                    energy_grid_(default_energy_grid),
                    beta_grid_(),
                    beta_cdf_grid_(default_cdf_grid),
                    alpha_cdf_grid_(default_cdf_grid) {
  // TODO
}

void Panglos::energy_grid(const std::vector<double>& egrid) {
  energy_grid_ = egrid;
}

const std::vector<double>& Panglos::energy_grid() const {
  return energy_grid_;
}

void Panglos::beta_grid(const std::vector<double>& bgrid) {
  beta_grid_ = bgrid;
}

const std::vector<double>& Panglos::beta_grid() const {
  return beta_grid_;
}

void Panglos::beta_cdf_grid(const std::vector<double>& bcgrid) {
  beta_cdf_grid_ = bcgrid;
}

const std::vector<double>& Panglos::beta_cdf_grid() const {
  return beta_cdf_grid_;
}

void Panglos::alpha_cdf_grid(const std::vector<double>& acgrid) {
  alpha_cdf_grid_ = acgrid;
}

const std::vector<double>& Panglos::alpha_cdf_grid() const {
  return alpha_cdf_grid_;
}

void Panglos::beta_max(double bmax) {
  beta_max_ = bmax;
}

double Panglos::beta_max() const {
  return beta_max_;
}

/* REFERENCES
 *
 * [1] A. T. Pavlou and W. Ji, “On-the-fly sampling of temperature-dependent
 *     thermal neutron scattering data for Monte Carlo simulations,” Ann Nucl
 *     Energy, vol. 71, pp. 411–426, 2014, doi: 10.1016/j.anucene.2014.04.028.  
 * 
 * */