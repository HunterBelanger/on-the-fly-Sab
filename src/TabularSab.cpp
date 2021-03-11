#include "TabularSab.hpp"
#include "constants.hpp"

using ScatteringFunction = section::Type<7,4>::TabulatedFunctions::ScatteringFunction;

TabularSab::TabularSab(section::Type<7,4>::TabulatedFunctions& TSL,
                       size_t indxT, double temp, double effTemp,
                       int lat, int lasym, int lln): symmetric_(lasym == 0),
                                                     temperature_(temp),
                                                     effTemperature_(effTemp),
                                                     SCTSab(temp, effTemp),
                                                     beta_(),
                                                     alpha_(),
                                                     data_(),
                                                     interpolations_() {
  //===========
  // Get the grid of beta values
  beta_ = TSL.betas();
  // Get the interpolations between beta values
  interpolations_ = TSL.interpolants();

  // Get all SaT records, one for each beta
  std::vector<ScatteringFunction> scatteringFuncs = TSL.S();

  // Make sure that beta_ and scatteringFuncs have the same size !
  if(beta_.size() != scatteringFuncs.size()) {
    throw std::runtime_error("Beta Grid and ScatteringFunctions of different sizes");
  }

  // The alpha grid is the same for all values of beta and T. We can then
  // get the alpha grid from the first scattering function.
  alpha_ = scatteringFuncs.front().alphas();

  // If lat = 1, then the alpha and beta grids need to be converted to
  // the true temperature, as they have been stored for room temp
  // T = 0.0253 eV.
  if(lat == 1) {
    double Troom = 0.0253 / kb;
    for(auto& beta : beta_) {
      beta *= (temperature_ / Troom);
    }
    for(auto& alpha : alpha_) {
      alpha *= (temperature_ / Troom);
    }
  }

  // Must now construct the Beta instance for each value of beta
  for(size_t b = 0; b < beta_.size(); b++) {
    std::vector<double> S = scatteringFuncs[b].S()[indxT];

    // Get the vector for the interpolations and boundaries along alpha, which is
    // the same for each beta, and temperatureA
    std::vector<long> alphaInterps = scatteringFuncs.front().interpolants();
    std::vector<long> alphaBounds = scatteringFuncs.front().boundaries();

    // Make sure that alpha_ and S have the same size !
    if(alpha_.size() != S.size()) {
      throw std::runtime_error("Alpha Grid and S grid of different sizes");
    }

    // If lln = 1, ln(S) is stored, and not S. We need to exponentiate
    if(lln == 1) {
      for(auto& s : S) s = std::exp(s);
    }

    // Apparently it is necessary to check for zeros in the Sa function, because
    // evaluators claim that the array is Log-Lin interporable, but we can't do
    // that when there are 0 values in the y array.
    fixData(alphaBounds, alphaInterps, S);

    // Create Tab1 for Sa
    Tab1 Sa = makeTab1({alphaBounds.begin(), alphaBounds.end()}, 
                       {alphaInterps.begin(), alphaInterps.end()},
                       {alpha_.begin(), alpha_.end()},
                       {S.begin(), S.end()});

    // Create beta instance
    data_.emplace_back(Beta(Sa, beta_[b]));
  }
}

void TabularSab::fixData(std::vector<long>& bounds, std::vector<long>& interps,
                         std::vector<double>& y) {

  // TODO
}