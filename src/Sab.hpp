#ifndef SAB_H
#define SAB_H

#include "constants.hpp"

#include <cmath>

class Sab {
public:
  Sab() {}
  virtual ~Sab() = default;

  virtual double operator()(double a, double b) const = 0;
  virtual double integrateAlpha(double aLow, double aHi, double b) const = 0;
  virtual double integrateAlphaExpBeta(double aLow, double aHi, double bLow, double bHi) const = 0;

  /**
   * @param E Incident energy in eV.
   * @param T Temperature in k.
   * @return Minimum possible value of beta for a given incident energy
   *         and temperature.
   */
  static double minBeta(double E, double T) {
    return -E / (kb * T);
  }

  /**
   * @param E Incident energy in eV.
   * @param b Value of beta.
   * @param T Temperature in k.
   * @param A Atomic Weight Ratio of nuclide.
   * @return Minimum possible value of alpha for a given incident energy, beta,
   *         temperature, and nuclide mass.
   */
  static double minAlpha(double E, double b, double T, double A) {
    return std::pow(std::sqrt(E) - std::sqrt(E+b*kb*T), 2.) / (A*kb*T);
  }

  /**
   * @param E Incident energy in eV.
   * @param b Value of beta.
   * @param T Temperature in k.
   * @param A Atomic Weight Ratio of nuclide.
   * @return Maximum possible value of alpha for a given incident energy, beta,
   *         temperature, and nuclide mass.
   */
  static double maxAlpha(double E, double b, double T, double A) {
    return std::pow(std::sqrt(E) + std::sqrt(E+b*kb*T), 2.) / (A*kb*T);
  }
  
};

#endif