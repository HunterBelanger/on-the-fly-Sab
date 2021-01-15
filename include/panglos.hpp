#ifndef PANGLOS_H
#define PANGLOS_H

#include <vector>

class Panglos {
  public:
    Panglos();
    ~Panglos() = default;

    //==========================================================================
    // No copy, move, or assignment constructors !
    Panglos(const Panglos&) = delete;
    Panglos(Panglos&&) = delete;
    Panglos& operator=(const Panglos&) = delete;
    Panglos& operator=(Panglos&&) = delete;

    //==========================================================================
    
    //==========================================================================
    // Getter and Setter Method for Sortage Grids
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
    double beta_max_; // Max energy transfer value to use (20 is used by default in Andrew's paper)

    // Grids along which Sab, PDF, and CDF tables are stored
    std::vector<double> temperature_; // [k]

    // Grids for storage of output
    std::vector<double> energy_grid_;
    std::vector<double> beta_grid_;
    std::vector<double> beta_cdf_grid_;
    std::vector<double> alpha_cdf_grid_;
};

#endif