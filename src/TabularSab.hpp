#ifndef TABULAR_SAB_H
#define TABULAR_SAB_H

#include <ENDFtk/section/7.hpp>

#include <vector>

#include "Sab.hpp"
#include "ShortCollisionTimeSab.hpp"
#include "Tab1.hpp"

class TabularSab : public Sab {
 public:
  TabularSab(section::Type<7, 4>::TabulatedFunctions& TSL, size_t indxT,
             double temp, double effTemp, int lat, int lasym, int lln);
  ~TabularSab() = default;

  double operator()(double a, double b) const override final {
    if (symmetric_ && b < 0.) b *= -1.;

    // Check that alpha and beta are in the given bounds.
    // If they aren't, use the short-collision-time approximation.
    if (a < alpha_.front() || a > alpha_.back() || b < beta_.front() ||
        b > beta_.back()) {
      return SCTSab(a, b);
    }

    // Get the right beta
    auto betaIt = std::lower_bound(beta_.begin(), beta_.end(), b);
    size_t betaIndex = std::distance(beta_.begin(), betaIt);
    if (*betaIt == b) return data_[betaIndex](a);

    double betaLow = beta_[betaIndex - 1];
    double betaHi = beta_[betaIndex];
    double SaBetaLow = data_[betaIndex - 1](a);
    double SaBetaHi = data_[betaIndex](a);

    // Get the interpolation rule for this interval in beta
    long interp = interpolations_[betaIndex - 1];

    if (interp == 1) {
      return ::Histogram::apply(b, betaLow, betaHi, SaBetaLow, SaBetaHi);
    } else if (interp == 2) {
      return ::LinearLinear::apply(b, betaLow, betaHi, SaBetaLow, SaBetaHi);
    } else if (interp == 3) {
      return ::LinearLogarithmic::apply(b, betaLow, betaHi, SaBetaLow,
                                        SaBetaHi);
    } else if (interp == 4) {
      return ::LogarithmicLinear::apply(b, betaLow, betaHi, SaBetaLow,
                                        SaBetaHi);
    } else if (interp == 5) {
      return ::LogarithmicLogarithmic::apply(b, betaLow, betaHi, SaBetaLow,
                                             SaBetaHi);
    } else {
      throw std::runtime_error(
          "Invalid interpolation along beta in TabularSab");
    }

    // Never gets here
    return 0.;
  }

  double integrateAlpha(double aLow, double aHi,
                        double b) const override final {
    if (symmetric_ && b < 0.) b *= -1.;

    if (b < beta_.front() || b > beta_.back()) {
      return SCTSab.integrateAlpha(aLow, aHi, b);
    } else {
      double integral = 0.;
      bool flipped = false;
      if (aLow > aHi) {
        flipped = true;
        std::swap(aLow, aHi);
      }

      // We are completely outsize of the alpha domain. Use the SCT
      if ((aLow < alpha_.front() && aHi < alpha_.front()) ||
          (aLow > alpha_.back() && aHi > alpha_.back())) {
        integral = SCTSab.integrateAlpha(aLow, aHi, b);
        if (flipped) {
          integral *= -1.;
        }
        return integral;
      }

      // Lower portion is outside of range
      if (aLow < alpha_.front()) {
        integral += SCTSab.integrateAlpha(aLow, alpha_.front(), b);
        aLow = alpha_.front();
      }

      // Upper portion is outside of range
      if (aHi > alpha_.back()) {
        integral += SCTSab.integrateAlpha(alpha_.back(), aHi, b);
        aHi = alpha_.back();
      }

      // Middle portion
      // Get the right beta
      auto betaIt = std::lower_bound(beta_.begin(), beta_.end(), b);
      size_t betaIndex = std::distance(beta_.begin(), betaIt);
      if (*betaIt == b) {
        integral += data_[betaIndex].integrate(aLow, aHi);
      } else {
        double betaLow = beta_[betaIndex - 1];
        double betaHi = beta_[betaIndex];
        double SaBetaLow = data_[betaIndex - 1].integrate(aLow, aHi);
        double SaBetaHi = data_[betaIndex].integrate(aLow, aHi);

        // Get the interpolation rule for this interval in beta
        long interp = interpolations_[betaIndex - 1];

        if (interp == 1) {
          integral +=
              ::Histogram::apply(b, betaLow, betaHi, SaBetaLow, SaBetaHi);
        } else if (interp == 2) {
          integral +=
              ::LinearLinear::apply(b, betaLow, betaHi, SaBetaLow, SaBetaHi);
        } else if (interp == 3) {
          integral += ::LinearLogarithmic::apply(b, betaLow, betaHi, SaBetaLow,
                                                 SaBetaHi);
        } else if (interp == 4) {
          integral += ::LogarithmicLinear::apply(b, betaLow, betaHi, SaBetaLow,
                                                 SaBetaHi);
        } else if (interp == 5) {
          integral += ::LogarithmicLogarithmic::apply(b, betaLow, betaHi,
                                                      SaBetaLow, SaBetaHi);
        } else {
          throw std::runtime_error(
              "Invalid interpolation along beta in TabularSab");
        }
      }

      // Once we have all portions, flip the sign if we flipped
      // the integration order.
      if (flipped) {
        integral *= -1.;
      }

      return integral;
    }

    // Never gets here
    return 0.;
  }

  double integrateAlphaExpBeta(double aLow, double aHi, double bLow,
                               double bHi) const override final {
    bool flipped = (bLow > bHi);
    if(flipped) {std::swap(bLow, bHi);}

    double integral = 0.;

    if ((bLow < beta_.front() && bHi < beta_.front()) ||
        (bLow > beta_.back() && bHi > beta_.back()) ) {
      integral += SCTSab.integrateAlphaExpBeta(aLow, aHi, bLow, bHi);
      if(flipped) integral *= -1.;
      return integral;
    }

    // Lower portion is outside of range
    if (bLow < beta_.front()) {
      integral += SCTSab.integrateAlphaExpBeta(aLow, aHi, bLow, beta_.front());
      bLow = beta_.front();
    }

    // Upper portion is outside of range
    if (bHi > beta_.back()) {
      integral += SCTSab.integrateAlphaExpBeta(aLow, aHi, beta_.back(), bHi);
      bHi = beta_.back();
    }

    // Get the middle range of the integral which is over the grid
    
    // This is the lambda which will be evaluated for each beta point in the
    // quadrature durring integration.
    auto expIntSa_at_b = [aLow, aHi, this](double b) {
      return std::exp(-b/2.) * this->integrateAlpha(aLow, aHi, b);
    };

    // Get the quadrature. Using 61 points, which can perfectly integrate
    // polynomials of up to order 184.
    GaussKronrodQuadrature<61> GKQ61;

    // The integral is returned as a pair of doubles. The first value is the
    // estimate of the integral, while the second value is the upper limit of
    // the error of the integral.
    std::pair<double,double> middle_integral = GKQ61.integrate(expIntSa_at_b, bLow, bHi);
    integral += middle_integral.first;


    if(flipped) {integral *= -1.;}
    return integral;
  }

  bool symmetric() const { return symmetric_; }

  std::vector<double> alphaGrid() const {return alpha_;}

  std::vector<double> betaGrid() const {return beta_;}

 private:
  bool symmetric_;
  double temperature_;
  double effTemperature_;
  ShortCollisionTimeSab SCTSab;
  std::vector<double> beta_;
  std::vector<double> alpha_;
  std::vector<Tab1> data_;
  std::vector<long> interpolations_;

  void fixData(std::vector<long>& bounds, std::vector<long>& interps,
               std::vector<double>& y);
};

#endif
