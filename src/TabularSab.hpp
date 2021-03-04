#ifndef TABULAR_SAB_H
#define TABULAR_SAB_H

#include <vector>

#include <ENDFtk/section/7.hpp>

#include "Sab.hpp"
#include "ShortCollisionTimeSab.hpp"
#include "Tab1.hpp"

class TabularSab : public Sab {
private:
  class Beta {
  public:
    Beta(Tab1 Sa, double b): beta_(b), Sa_(Sa) {}
    ~Beta() = default;

    double operator()(double a) const {
      return this->Sa_(a);
    }
    double integrate(double aLow, double aHi) const {
      return this->Sa_.integrate(aLow, aHi);
    }
    double value() const {return this->beta_;}

  private:
    double beta_;
    Tab1 Sa_; // Scattering law as a function of alpha
  };

public:
  TabularSab(section::Type<7,4>::TabulatedFunctions& TSL, size_t indxT, double temp, double effTemp, int lat, int lasym, int lln);
  ~TabularSab() = default;

  double operator()(double a, double b) const override final;
  double integrateAlpha(double aLow, double aHi, double b) const override final;
  double integrateAlphaExpBeta(double aLow, double aHi, double bLow, double bHi) const override final;

  bool symmetric() const;

private:
  bool symmetric_;
  double temperature_;
  double effTemperature_;
  ShortCollisionTimeSab SCTSab;
  std::vector<double> beta_;
  std::vector<double> alpha_;
  std::vector<Beta> data_;
  // TODO interpolation along beta
};

#endif
