#include "TabularSab.hpp"

TabularSab::TabularSab(section::Type<7,4>::TabulatedFunctions& TSL,
                       size_t indxT, double temp, double effTemp,
                       int lat, int lasym, int lln): symmetric_(lasym == 0),
                                                     temperature_(temp),
                                                     effTemperature_(effTemp),
                                                     SCTSab(temp, effTemp),
                                                     beta_(),
                                                     alpha_(),
                                                     data_() {
 //===========
  beta_ = TSL.betas();
  
}
