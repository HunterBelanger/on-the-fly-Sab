#ifndef SAB_H
#define SAB_H

class Sab {
public:
  Sab(double AWR, double temp): AWR_(AWR), temperature_(temp) {}
  virtual ~Sab() = default;

  virtual double operator()(double a, double b) const = 0;
  virtual double integrateAlpha(double aLow, double aHi, double b) const = 0;
  
  double temperature() const {
    return this->temperature_;
  }

  // Functions to calculate the integration limits for the
  // generation of PDFs.
  double minEnergyTransfer(double E) const; // beta_min
  double minMomentumTransfer(double E, double b) const; // alpha_min
  double maxMomentumTransfer(double E, double b) const; // alpha_max

protected:
  double AWR_;
  double temperature_;
};

#endif