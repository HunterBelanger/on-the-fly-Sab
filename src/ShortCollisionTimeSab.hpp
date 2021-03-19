#ifndef SHORT_COLLISION_TIME_SAB_H
#define SHORT_COLLISION_TIME_SAB_H

#include "Sab.hpp"
#include "GaussKronrod.hpp"
#include "constants.hpp"

class ShortCollisionTimeSab : public Sab {
 public:
  ShortCollisionTimeSab(double T, double Teff) : R(Teff / T) {}
  ~ShortCollisionTimeSab() = default;

  double operator()(double a, double b) const override final {
    b = std::abs(b);
    double numerator = std::exp(-((a - b) * (a - b) / (4. * a * R)) - 0.5 * b);
    double denominator = std::sqrt(4 * pi * a * R);
    return numerator / denominator;
  }

  double integrateAlpha(double aLow, double aHi,
                        double b) const override final {
    return this->indefiniteIntegralAlpha(aHi, b) -
           this->indefiniteIntegralAlpha(aLow, b);
  }

  double integrateAlphaExpBeta(double aLow, double aHi, double bLow,
                               double bHi) const override final {

    // Integration along beta is performed using the Gauss-Kronrod quadrature
    // while the alpha integration is still done analytically.

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
    std::pair<double,double> integral = GKQ61.integrate(expIntSa_at_b, bLow, bHi);

    return integral.first;
  }

 private:
  double R;

  double indefiniteIntegralAlpha(double a, double b) const {
    b = std::abs(b);
    double C1 = 0.5 * std::exp(-0.5 * b);
    double T1 = std::erf((a - b) / (2. * std::sqrt(a * R)));
    double C2 = std::exp(b / R);
    double T2 = std::erf((a + b) / (2. * std::sqrt(a * R)));
    return C1 * (T1 + C2 * (T2 - 1.) + 1.);
  }
};

#endif