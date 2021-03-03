#include <Panglos.hpp>

#include <ENDFtk/details.hpp>

// Private header files
#include "Sab.hpp"
#include "Tabular2D.hpp"

extern std::vector<double> defaultEnergyGrid;
extern std::vector<double> defaultCDFGrid;
extern std::vector<double> defaultTempertures;

Panglos::Panglos(file::Type<7>& mf7): mf7(mf7),
                                      LAT(),
                                      LASYM(),
                                      LLN(),
                                      AWR(),
                                      maxBeta_(20.),
                                      temps_(),
                                      TSLs_(),
                                      alphaCDFs_(),
                                      betaCDFs_(),
                                      energyGrid_(defaultEnergyGrid),
                                      betaGrid_(),
                                      betaCDFGrid_(defaultCDFGrid),
                                      alphaCDFGrid_(defaultCDFGrid) {

  // No idea what the '_c' is for, but it doesn't work without it
  section::Type<7, 4> mt4 = this->mf7.section(4_c);
  this->AWR = mt4.AWR();
  this->LAT = mt4.LAT();
  this->LASYM = mt4.LASYM();

  section::Type<7,4>::ScatteringLawConstants slConstants = mt4.constants();
  this->LLN = slConstants.LLN();

  // Check to see if using an analytic function or not
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

//==============================================================================
// Private Initialization Methods

void Panglos::initializeScatteringLaws() {

}

void Panglos::initializeAnalyticScatteringLaws() {
  // Since we are using analytic laws, no temperature grid is provided. As such
  // we use the default temperature grid.
  temps_ = defaultTempertures;

}

void Panglos::initializeTabularScatteringLaws() {
  // If LLN is set, then ln(S) is stored, and not S ! When this is the case,
  // must exponentiate all of the values so that integrals are calculated
  // correctly. In addition, one must also change the interpolation rule !!
  // If LogLog was used for (a,ln(S)), then this becomes LogLin once
  // exponentiated for (a,S) !!
  // LogLog -> LogLin
  // LinLog -> LinLin
  // I don't know what to do with other interpolation rules yet...
  // Maybe I will just leave them unchanged, so linear in S ?
}