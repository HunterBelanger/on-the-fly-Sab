#ifndef FREE_GAS_SAB_H
#define FREE_GAS_SAB_H

#include "Sab.hpp"
#include "constants.hpp"

#include <cmath>

class FreeGasSab : public Sab {
public:
  FreeGasSab() {}
  ~FreeGasSab() = default;

  double operator()(double a, double b) const override final {
    return (1. / std::sqrt(4.*pi*a))*std::exp(-(a*a + b*b)/(4. * a));
  }

  double integrateAlpha(double aLow, double aHi, double b) const override final {
    return this->indefiniteIntegralAlpha(aHi,b) - this->indefiniteIntegralAlpha(aLow,b);
  }

  double integrateAlphaExpBeta(double aLow, double aHi, double bLow, double bHi) const override final {
    return 0.5*(this->indefiniteIntegralAlphaExpBeta1(aLow, aHi, bHi) -
                this->indefiniteIntegralAlphaExpBeta1(aLow, aHi, bLow));
  }

private:
  double indefiniteIntegralAlpha(double a, double b) const {
    b = std::abs(b);
    double Coeff = 0.5 * std::exp(-0.5*b);
    double T1 = -std::erf((0.5*b/std::sqrt(a))-(0.5*std::sqrt(a)));
    double T2 = std::exp(b)*(std::erf((b + a)/(2.*std::sqrt(a))) - 1.);
    return Coeff * (T1 + T2 + 1.);
  }

  double indefiniteIntegralAlphaExpBeta1(double aLow, double aHi, double b) const {
    return this->indefiniteIntegralAlphaExpBeta2(aHi,b) -
           this->indefiniteIntegralAlphaExpBeta2(aLow, b);
  }

  double indefiniteIntegralAlphaExpBeta2(double a, double b) const {
    double abs_b = std::abs(b);
    double C1 = 0.5 * std::exp(0.5*(-abs_b - b));
    double T1 = (b - abs_b - 2.)*std::erf((a - abs_b)/(2.*std::sqrt(a)));
    double T2 = (abs_b + b - 2.)*(std::exp(abs_b)*std::erf((a + abs_b)/(2.*std::sqrt(a)))-std::exp(abs_b) + 1.);
    double C2 = a * std::erf((a+b)/(2.*std::sqrt(a)));
    double C3 = (2.*std::sqrt(a)*std::exp(-(a+b)*(a+b)/(4.*a)))/(std::sqrt(pi));
    return C1*(T1+T2) + C2 + C3;
  }

};

#endif