#include <panglos.hpp>

// Private header files
#include "sab.hpp"
#include "tabular_2d.hpp"

extern std::vector<double> default_energy_grid;
extern std::vector<double> default_cdf_grid;
 

Panglos::Panglos(ENDFtk::file::Type<7>& mf7): mf7(mf7),
                                              beta_max_(20.),
                                              temps_(),
                                              tsls_(),
                                              alpha_pdfs_(),
                                              alpha_cdfs_(),
                                              beta_pdfs_(),
                                              beta_cdfs_(),
                                              energy_grid_(default_energy_grid),
                                              beta_grid_(),
                                              beta_cdf_grid_(default_cdf_grid),
                                              alpha_cdf_grid_(default_cdf_grid) {
  // TODO
}

//==============================================================================
// Getter and Setter Methods
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