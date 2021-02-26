#ifndef PANGLOS_H
#define PANGLOS_H

#include <vector>
#include <memory>

// Include ONLY this file !! Otherwise compilation will take forever !
#include <ENDFtk/file/7.hpp>
using namespace njoy;

class Sab;
class Tabular2D;

class Panglos {
public:
  Panglos(ENDFtk::file::Type<7>& mf7);
  ~Panglos() = default;

  //==========================================================================
  // No copy, move, or assignment constructors !
  Panglos(const Panglos&) = delete;
  Panglos(Panglos&&) = delete;
  Panglos& operator=(const Panglos&) = delete;
  Panglos& operator=(Panglos&&) = delete;

  //==========================================================================
  
  //==========================================================================
  // Getter and Setter Methods for Sortage Grids
  void energyGrid(const std::vector<double>& egrid);
  const std::vector<double>& energyGrid() const;
  
  void betaGrid(const std::vector<double>& bgrid);
  const std::vector<double>& betaGrid() const;

  void betaCDFGrid(const std::vector<double>& bcgrid);
  const std::vector<double>& betaCDFGrid() const;
  
  void alphaCDFGrid(const std::vector<double>& acgrid);
  const std::vector<double>& alphaCDFGrid() const;

  double maxEnergyTransfer() const;
  void maxEnergyTransfer(double bmax);

private:
  ENDFtk::file::Type<7> mf7;
  
  // Max energy transfer value to use (20 is used by default in Andrew's paper)
  double maxEnergyTransfer_;

  // Grids along which Sab, PDF, and CDF tables are stored
  std::vector<double> temps_; // [k]
  std::vector<std::unique_ptr<Sab>> TSLs_;
  std::vector<std::unique_ptr<Tabular2D>> alphaPDFs_;
  std::vector<std::unique_ptr<Tabular2D>> alphaCDFs_;
  std::vector<std::unique_ptr<Tabular2D>> betaPDFs_;
  std::vector<std::unique_ptr<Tabular2D>> betaCDFs_;

  // Grids for storage of output
  std::vector<double> energyGrid_;
  std::vector<double> betaGrid_;
  std::vector<double> betaCDFGrid_;
  std::vector<double> alphaCDFGrid_;
};

#endif