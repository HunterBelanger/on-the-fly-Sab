#ifndef TABULAR_SAB_H
#define TABULAR_SAB_H

#include <vector>

#include "Sab.hpp"
#include "Tab1.hpp"

class TabularSab : public Sab {
  private:
    class Beta {
      public:
        Beta();
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
    TabularSab(); // This will have to take the file, and the temperature index
    ~TabularSab() = default;

    double operator()(double a, double b) const override final;
    double integrateAlpha(double aLow, double aHi, double b) const override final;

    bool symmetric() const;

  private:
    bool symmetric_;
    double effTemperature_;
    std::vector<double> beta_;
    std::vector<double> alpha_;
    std::vector<Beta> data_;
    // TODO interpolation along beta
};

#endif