#ifndef SHORT_COLLISION_TIME_SAB_H
#define SHORT_COLLISION_TIME_SAB_H

#include "Sab.hpp"
#include "constants.hpp"

class ShortCollisionTimeSab : public Sab {
public:
  ShortCollisionTimeSab(double T, double Teff): R(Teff / T) {}
  ~ShortCollisionTimeSab() = default;

  double operator()(double a, double b) const override final {
    b = std::abs(b);
    double numerator = std::exp(-((a-b)*(a-b)/(4.*a*R)) - 0.5*b);
    double denominator = std::sqrt(4*pi*a*R);
    return numerator / denominator;
  }

  double integrateAlpha(double aLow, double aHi, double b) const override final {
    return this->indefiniteIntegralAlpha(aHi, b) - this->indefiniteIntegralAlpha(aLow, b);
  }

  double integrateAlphaExpBeta(double aLow, double aHi, double bLow, double bHi) const override final;

private:
  double R;

  double indefiniteIntegralAlpha(double a, double b) const {
    b = std::abs(b);
    double C1 = 0.5*std::exp(-0.5*b);
    double T1 = std::erf((a - b)/(2.*std::sqrt(a*R)));
    double C2 = std::exp(b/R);
    double T2 = std::erf((a + b)/(2.*std::sqrt(a*R)));
    return C1 * (T1 + C2 * (T2 - 1.) + 1.);
  }

};

#endif