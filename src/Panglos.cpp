#include <Panglos.hpp>

// Private header files
#include "Sab.hpp"
#include "Tabular2D.hpp"

extern std::vector<double> defaultEnergyGrid;
extern std::vector<double> defaultCDFGrid;
 

Panglos::Panglos(ENDFtk::file::Type<7>& mf7): mf7(mf7),
                                              maxBeta_(20.),
                                              temps_(),
                                              TSLs_(),
                                              alphaPDFs_(),
                                              alphaCDFs_(),
                                              betaPDFs_(),
                                              betaCDFs_(),
                                              energyGrid_(defaultEnergyGrid),
                                              betaGrid_(),
                                              betaCDFGrid_(defaultCDFGrid),
                                              alphaCDFGrid_(defaultCDFGrid) {
  // TODO
}

//==============================================================================
// Getter and Setter Methods
void Panglos::energyGrid(const std::vector<double>& egrid) {
  this->energyGrid_ = egrid;
}

const std::vector<double>& Panglos::energyGrid() const {
  return this->energyGrid_;
}

void Panglos::betaGrid(const std::vector<double>& bgrid) {
  this->betaGrid_ = bgrid;
}

const std::vector<double>& Panglos::betaGrid() const {
  return this->betaGrid_;
}

void Panglos::betaCDFGrid(const std::vector<double>& bcgrid) {
  this->betaCDFGrid_ = bcgrid;
}

const std::vector<double>& Panglos::betaCDFGrid() const {
  return this->betaCDFGrid_;
}

void Panglos::alphaCDFGrid(const std::vector<double>& acgrid) {
  this->alphaCDFGrid_ = acgrid;
}

const std::vector<double>& Panglos::alphaCDFGrid() const {
  return this->alphaCDFGrid_;
}

void Panglos::maxBeta(double bmax) {
  this->maxBeta_ = bmax;
}

double Panglos::maxBeta() const {
  return this->maxBeta_;
}