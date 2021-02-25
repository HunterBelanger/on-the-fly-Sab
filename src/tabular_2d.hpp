#ifndef TABULAR_2D_H
#define TABULAR_2D_H

#include <vector>

// Include logging library so that errors or warnings can be written,
// should integration or evaluation be attempted for invalid values.
class Tabular2D {
  public:
    Tabular2D(); // Will take x and y in, and default value for all data elements
    Tabular2D(size_t nx, size_t ny);
    // On construction, should assure that both x and y are sorted !

    ~Tabular2D() = default;

    double& operator()(size_t ix, size_t iy);
    const double& operator()(size_t ix, size_t iy) const;

    double evaluate(double x, double y) const;
    double integrate(double xlow, double xhi, double ylow, double yhi) const;

    double& x(size_t i) { return x_[i]; }
    const double& x(size_t i) const { return x_[i]; }

    double& y(size_t i) { return y_[i]; }
    const double& y(size_t i) const { return y_[i]; }

    const std::vector<double>& x_grid() const { return x_; }
    const std::vector<double>& y_grid() const { return y_; }

    size_t x_size() const { return x_.size(); }
    size_t y_size() const { return y_.size(); }

    double x_min() const { return x_.front(); }
    double x_max() const { return x_.back(); }
    double y_min() const { return y_.front(); }
    double y_max() const { return y_.back(); }

  private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> data_; 
};

#endif