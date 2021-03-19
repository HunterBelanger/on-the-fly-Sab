#include "TabularSab.hpp"

#include "constants.hpp"

using ScatteringFunction =
    section::Type<7, 4>::TabulatedFunctions::ScatteringFunction;

TabularSab::TabularSab(section::Type<7, 4>::TabulatedFunctions& TSL,
                       size_t indxT, double temp, double effTemp, int lat,
                       int lasym, int lln)
    : symmetric_(lasym == 0),
      temperature_(temp),
      effTemperature_(effTemp),
      SCTSab(temp, effTemp),
      beta_(),
      alpha_(),
      data_(),
      interpolations_() {
  //============================================================================
  // Get the grid of beta values
  beta_ = TSL.betas();
  // Get the interpolations between beta values
  interpolations_ = TSL.interpolants();

  // Get all SaT records, one for each beta
  std::vector<ScatteringFunction> scatteringFuncs = TSL.S();

  // Make sure that beta_ and scatteringFuncs have the same size !
  if (beta_.size() != scatteringFuncs.size()) {
    throw std::runtime_error(
        "Beta Grid and ScatteringFunctions of different sizes");
  }

  // The alpha grid is the same for all values of beta and T. We can then
  // get the alpha grid from the first scattering function.
  alpha_ = scatteringFuncs.front().alphas();

  // If lat = 1, then the alpha and beta grids need to be converted to
  // the true temperature, as they have been stored for room temp
  // T = 0.0253 eV.
  if (lat == 1) {
    double Troom = 0.0253 / kb;
    for (auto& beta : beta_) {
      beta *= (temperature_ / Troom);
    }
    for (auto& alpha : alpha_) {
      alpha *= (temperature_ / Troom);
    }
  }

  // Get the vector for the interpolations and boundaries along alpha, which is
  // the same for each beta, and temperatureA
  std::vector<long> alphaInterpsCommon = scatteringFuncs.front().interpolants();
  std::vector<long> alphaBoundsCommon = scatteringFuncs.front().boundaries();

  // If ln(S) is stored, we need to change the interpolation rules.
  // ATTENTION ! The way the ENDF manual describes this is wrong !
  // Written bellow are the proper transformations.
  // LinLog(3) -> LogLog(5)
  // LinLin(2) -> LogLin(4)
  // I don't know what to do with other interpolation rules yet...
  // Maybe I will just leave them unchanged, so linear in S ?
  // In theory, I don't think an evaluator would be allowed to use
  // one of the other interpolation rules when storing ln(S), as they
  // would have the same issue. I should probably add a check here ?
  if (lln == 1) {
    for (auto& interp : alphaInterpsCommon) {
      if (interp == 3) {
        interp = 5;
      } else if (interp == 2) {
        interp = 4;
      } else {
        throw std::runtime_error("Cannot have LLN=1 with given interp law");
      }
    }
  }

  // Must now construct the Beta instance for each value of beta
  for (size_t b = 0; b < beta_.size(); b++) {
    std::vector<double> S = scatteringFuncs[b].S()[indxT];

    // Need to make copies of these two vectors as they typically need
    // to be modified depending on the zeros in the fixData method.
    std::vector<long> alphaInterps = alphaInterpsCommon;
    std::vector<long> alphaBounds = alphaBoundsCommon;

    // Make sure that alpha_ and S have the same size !
    if (alpha_.size() != S.size()) {
      throw std::runtime_error("Alpha Grid and S grid of different sizes");
    }

    // If lln = 1, ln(S) is stored, and not S. We need to exponentiate.
    // The interpolation rules were already changed above.
    if (lln == 1) {
      for (auto& s : S) {
        s = std::exp(s);
      }
    }

    // Apparently it is necessary to check for zeros in the Sa function, because
    // evaluators claim that the array is Log-Lin interporable, but we can't do
    // that when there are 0 values in the y array.
    fixData(alphaBounds, alphaInterps, S);

    // Create Tab1 for Sa
    Tab1 Sa = makeTab1({alphaBounds.begin(), alphaBounds.end()},
                       {alphaInterps.begin(), alphaInterps.end()},
                       {alpha_.begin(), alpha_.end()}, {S.begin(), S.end()});

    // Save table to data
    data_.push_back(Sa);
  }
}

void TabularSab::fixData(std::vector<long>& bounds, std::vector<long>& interps,
                         std::vector<double>& y) {
  // Only need to fix data is one region is given, and the region is a Log-x
  // interpolation type (y is log)
  if ((bounds.size() == 1) && ((interps[0] == 4) || (interps[0] == 5))) {
    long origInterp = interps.front();

    bounds.clear();
    interps.clear();

    for (long i = 0; i < static_cast<long>(y.size()) - 1; i++) {
      if (y[i] == 0. && y[i + 1] != 0.) {
        if (i < static_cast<long>(y.size()) - 2 && y[i + 2] == 0.) {
          // Here we have a segment where
          // [ ..., 0,  X,    0, ...]
          //        i,  i+1,  i+2
          // In this case, since we only have a single non-zero value, we
          // keep the linear inerpolation and don't add a break point but
          // we do need to advance i by 1
          i += 1;
        } else {
          // Here we have a segment where
          // [ ...., 0, X, Y, ...]
          // In this case, we record the break point !
          bounds.push_back(i + 2);
          interps.push_back(2);  // Ending a lin-lin region
        }
      } else if (y[i] != 0. && y[i + 1] == 0.) {
        // Here we have a segment
        // [ ...., X, Y, 0, ...]
        bounds.push_back(i + 1);
        interps.push_back(origInterp);
      }
    }

    // Add the end interval info, which isn't treated by the
    // for-loop above.
    bounds.push_back(static_cast<long>(y.size()));
    if (y[y.size() - 2] == 0. || y[y.size() - 1] == 0.) {
      interps.push_back(2);
    } else {
      interps.push_back(origInterp);
    }
  }

  if (!std::is_sorted(bounds.begin(), bounds.end())) {
    throw std::runtime_error("Bad fixing of data");
  }
}