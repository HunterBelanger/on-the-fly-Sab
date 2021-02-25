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
    void energy_grid(const std::vector<double>& egrid);
    const std::vector<double>& energy_grid() const;
    
    void beta_grid(const std::vector<double>& bgrid);
    const std::vector<double>& beta_grid() const;

    void beta_cdf_grid(const std::vector<double>& bcgrid);
    const std::vector<double>& beta_cdf_grid() const;
    
    void alpha_cdf_grid(const std::vector<double>& acgrid);
    const std::vector<double>& alpha_cdf_grid() const;

    double beta_max() const;
    void beta_max(double bmax);

  private:
    ENDFtk::file::Type<7> mf7;
    
    // Max energy transfer value to use (20 is used by default in Andrew's paper)
    double beta_max_;

    // Grids along which Sab, PDF, and CDF tables are stored
    std::vector<double> temps_; // [k]
    std::vector<std::unique_ptr<Sab>> tsls_;
    std::vector<std::unique_ptr<Tabular2D>> alpha_pdfs_;
    std::vector<std::unique_ptr<Tabular2D>> alpha_cdfs_;
    std::vector<std::unique_ptr<Tabular2D>> beta_pdfs_;
    std::vector<std::unique_ptr<Tabular2D>> beta_cdfs_;

    // Grids for storage of output
    std::vector<double> energy_grid_;
    std::vector<double> beta_grid_;
    std::vector<double> beta_cdf_grid_;
    std::vector<double> alpha_cdf_grid_;
};

#endif