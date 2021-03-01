#ifndef GUASSKRONROD_H
#define GAUSSKRONROD_H

#include <cmath>
#include <functional>
#include <vector>

template<size_t NL, size_t NK>
struct GaussKronrodQuadrature {
  static std::pair<double, double> integrate(std::function<double(double)> func, double xLow, double xHi) {
    // First calculate integral using the Gauss-Legendre quadrature
    double glIntegral = 0.;
    double gkIntegral = 0;
    for(size_t i = 0; i < glWeights.size(); i++) {
      double func_x = 0.;

      // Get the positive point
      double xi = abscissae[i];
      double x = 0.5*(xHi - xLow)*xi + 0.5*(xLow + xHi);
      func_x = func(x);
      glIntegral += func_x * glWeights[i];
      gkIntegral += func_x * weights[i];

      // Get the negative point
      if(xi != 0.) {
        xi *= -1.;
        double x = 0.5*(xHi - xLow)*xi + 0.5*(xLow + xHi);
        func_x = func(x);
        glIntegral += func_x * glWeights[i];
        gkIntegral += func_x * weights[i];
      }
    }

    // Now calculate integral using the Gauss-Kronrod quadrature
    // Get contribution from Gauss-Kronrod points
    for(size_t i = glWeights.size(); i < weights.size(); i++) {
      // Get the positive point
      double xi = abscissae[i];
      double x = 0.5*(xHi - xLow)*xi + 0.5*(xLow + xHi);
      gkIntegral += func(x) * weights[i];

      // Get the negative point
      if(xi != 0.) {
        xi *= -1.;
        double x = 0.5*(xHi - xLow)*xi + 0.5*(xLow + xHi);
        gkIntegral += func(x) * weights[i];
      }
    }

    glIntegral *= 0.5*(xHi - xLow);
    gkIntegral *= 0.5*(xHi - xLow);

    // Calculate the error
    double err = std::abs(glIntegral - gkIntegral) / gkIntegral;

    return {gkIntegral, err};
  }

  static constexpr size_t order() {return 3*NK + 1;}

  static const std::vector<double> abscissae;
  static const std::vector<double> weights;
  static const std::vector<double> glWeights;
};

extern template struct GaussKronrodQuadrature< 7, 15>;
extern template struct GaussKronrodQuadrature<10, 21>;
extern template struct GaussKronrodQuadrature<15, 31>;
extern template struct GaussKronrodQuadrature<20, 41>;
extern template struct GaussKronrodQuadrature<25, 51>;

#endif