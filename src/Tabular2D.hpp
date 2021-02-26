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
    double integrate(double xLow, double xHi, double yLow, double yHi) const;

    double& x(size_t i) { return this->x_[i]; }
    const double& x(size_t i) const { return this->x_[i]; }

    double& y(size_t i) { return this->y_[i]; }
    const double& y(size_t i) const { return this->y_[i]; }

    const std::vector<double>& xGrid() const { return this->x_; }
    const std::vector<double>& yGrid() const { return this->y_; }

    size_t xSize() const { return this->x_.size(); }
    size_t ySize() const { return this->y_.size(); }

    double xMin() const { return this->x_.front(); }
    double xMax() const { return this->x_.back(); }
    double yMin() const { return this->y_.front(); }
    double yMax() const { return this->y_.back(); }

  private:
    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> data_; 
};

#endif