#ifndef SAB_H
#define SAB_H

#include <vector>

#include <interpolation.hpp>
using namespace njoy::interpolation;

using Law1 = Table<table::Type<Histogram,
                               table::search::Binary,
                               table::discontinuity::TakeLeft,
                               std::vector<double>,
                               std::vector<double>>,
                   table::left::interval::Throws,
                   table::right::interval::Throws>;

using Law2 = Table<table::Type<LinearLinear,
                               table::search::Binary,
                               table::discontinuity::TakeLeft,
                               std::vector<double>,
                               std::vector<double>>,
                   table::left::interval::Throws,
                   table::right::interval::Throws>;

using Law3 = Table<table::Type<LinearLogarithmic,
                               table::search::Binary,
                               table::discontinuity::TakeLeft,
                               std::vector<double>,
                               std::vector<double>>,
                   table::left::interval::Throws,
                   table::right::interval::Throws>;

using Law4 = Table<table::Type<LogarithmicLinear,
                               table::search::Binary,
                               table::discontinuity::TakeLeft,
                               std::vector<double>,
                               std::vector<double>>,
                   table::left::interval::Throws,
                   table::right::interval::Throws>;

using Law5 = Table<table::Type<LogarithmicLogarithmic,
                               table::search::Binary,
                               table::discontinuity::TakeLeft,
                               std::vector<double>,
                               std::vector<double>>,
                   table::left::interval::Throws,
                   table::right::interval::Throws>;

using ENDFvariant = Table<table::Variant<Law1, Law2, Law3, Law4, Law5>>;
using Tab1 = Table<table::Vector<ENDFvariant>>;

class Sab {
  public:

    class Beta {
      public:
        Beta();
        ~Beta() = default;

        double operator()(double a) const {
          return Sa_(a);
        }
        double integrate(double a_low, double a_hi) const {
          return Sa_.integrate(a_low, a_hi);
        }
        double value() const {return beta_;}

      private:
        double beta_;
        Tab1 Sa_; // Scattering law as a function of alpha
    };

  public:
    Sab(); // This will have to take the file, and the temperature index
    ~Sab() = default;

    double operator()(double a, double b) const;
    double integrate_alpha(double b, double a_low, double a_hi) const;

    const Beta& operator[](size_t i) const;

    const std::vector<double>& beta() const;
    double beta(size_t i) const;
    size_t nbeta() const;

    const std::vector<double>& alpha() const;
    double alpha(size_t i) const;
    size_t nalpha() const;

    double temperature() const;
    bool symmetric() const;

    // Functions to calculate the integration limits for the
    // generation of PDFs.
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