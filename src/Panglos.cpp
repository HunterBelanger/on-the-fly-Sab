#include <Panglos.hpp>
#include <boost/hana.hpp>  // Needed for the _c literal for constructing mt4
#include <iostream>

// Private header files
#include "Sab.hpp"
#include "Tab1.hpp"
#include "TabularSab.hpp"

extern std::vector<double> defaultEnergyGrid;
extern std::vector<double> defaultCDFGrid;
extern std::vector<double> defaultTempertures;

Panglos::Panglos(file::Type<7>& mf7)
    : mf7(mf7),
      mt4(mf7.section(4_c)),
      LAT(),
      LASYM(),
      LLN(),
      AWR(),
      maxBeta_(20.),
      temps_(),
      TSLs_(),
      // alphaCDFs_(),
      // betaCDFs_(),
      energyGrid_(defaultEnergyGrid),
      betaGrid_(),
      betaCDFGrid_(defaultCDFGrid),
      alphaCDFGrid_(defaultCDFGrid) {
  this->AWR = this->mt4.AWR();
  this->LAT = this->mt4.LAT();
  this->LASYM = this->mt4.LASYM();

  section::Type<7, 4>::ScatteringLawConstants slConstants =
      this->mt4.constants();
  this->LLN = slConstants.LLN();

  initializeScatteringLaws();
}

// Must put destructor here in cpp due to PIMPL Sab in header
Panglos::~Panglos() {}

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

const std::vector<double>& Panglos::betaGrid() const { return this->betaGrid_; }

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

void Panglos::maxBeta(double bmax) { this->maxBeta_ = bmax; }

double Panglos::maxBeta() const { return this->maxBeta_; }

//==============================================================================
// Private Initialization Methods

void Panglos::initializeScatteringLaws() {
  // Get the actual scattering law, which is provided as a variant of the
  // analytic or tabular types
  auto scatteringLaw = this->mt4.scatteringLaw();

  // See if the variant holds a tabulated or analytic Sab
  if (std::holds_alternative<section::Type<7, 4>::TabulatedFunctions>(
          scatteringLaw)) {
    initializeTabularScatteringLaws(
        std::get<section::Type<7, 4>::TabulatedFunctions>(scatteringLaw));
  } else {
    // Variant holds an analytic law for primary nuclide
    initializeAnalyticScatteringLaws(
        std::get<section::Type<7, 4>::AnalyticalFunctions>(scatteringLaw));
  }
}

void Panglos::initializeAnalyticScatteringLaws(
    section::Type<7, 4>::AnalyticalFunctions& /*tsl*/) {
  std::cout
      << " HELP ME ! I don't know what to do with analytic functions yet !\n";
}

void Panglos::initializeTabularScatteringLaws(
    section::Type<7, 4>::TabulatedFunctions& tsl) {
  // Get list of all provided temperatures from the scattering law
  auto lawsForAllTemps = tsl.S();
  temps_ = lawsForAllTemps[0].T();

  section::Type<7, 4>::EffectiveTemperature rawEffectiveTemp =
      mt4.principalEffectiveTemperature();

  // Get Tab1 for the effective temperature
  Tab1 effectiveTemp =
      makeTab1(rawEffectiveTemp.boundaries(), rawEffectiveTemp.interpolants(),
               rawEffectiveTemp.TMOD(), rawEffectiveTemp.TEFF());

  // Go through all temps and build each scatteriing law
  size_t Tindx = 0;
  for (const auto& temp : temps_) {
    std::cout << " Mod Temp = " << temp
              << ", Eff Temp = " << effectiveTemp(temp) << std::endl;

    double effTemp = effectiveTemp(temp);

    TSLs_.emplace_back(std::make_unique<TabularSab>(tsl, Tindx, temp, effTemp,
                                                    LAT, LASYM, LLN));
    Tindx++;
  }
}
