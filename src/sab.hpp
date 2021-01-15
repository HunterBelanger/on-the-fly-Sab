#ifndef SAB_H
#define SAB_H

#include <vector>

class Sab {
  public:

    class Beta {
      public:
        Beta();
        ~Beta() = default;

        double operator()(double a) const;
        double integrate(double a_low, double a_hi) const;

      private:
        double beta_;
        // TODO make Tab1 structure from interpolation library
    };

  public:
    Sab(); // This will have to take the file, and the temperature index
    ~Sab() = default;

    double operator()(double a, double b) const;
    double integrate(double a_low, double a_hi, double b_low, double b_hi) const;

    const Beta& operator[](size_t i) const;

    const std::vector<double>& beta() const;
    double beta(size_t i) const;
    size_t nbeta() const;

    const std::vector<double>& alpha() const;
    double alpha(size_t i) const;
    size_t nalpha() const;

    double temperature() const;
    bool symmetric() const;

    double alpha_min(double E, double b) const;
    double alpha_max(double E, double b) const;
    double beta_min(double E) const;

  private:
    double A_;
    double temperature_;
    bool symmetric_;
    std::vector<double> beta_;
    std::vector<double> alpha_;
    std::vector<Beta> data_;
    // TODO interpolation along beta
};

#endif