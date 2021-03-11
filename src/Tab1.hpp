#ifndef TAB1_H
#define TAB1_H

#include <interpolation.hpp>
using namespace njoy::interpolation;

#include <ENDFtk/types.hpp>
using namespace njoy::ENDFtk;

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

// Function to generate a Tab1 from the breakpoints, interpolations,
// x array, and y array
inline Tab1 makeTab1(LongRange breakpoints, LongRange interpolations, DoubleRange x, DoubleRange y) {
  std::vector<ENDFvariant> core;

  auto brkSize = ranges::size(breakpoints);
  size_t low = 0;
  size_t hi = 0;
  for(size_t i = 0; i < brkSize; i++) {
    int law = *(interpolations.begin()+i);
    hi = *(breakpoints.begin()+i);

    std::vector<double> xInterval = {x.begin()+low, x.begin()+hi};
    std::vector<double> yInterval = {y.begin()+low, y.begin()+hi};

    if(law == 1) {
      core.push_back(Law1(std::move(xInterval), std::move(yInterval)));
    } else if(law == 2) {
      core.push_back(Law2(std::move(xInterval), std::move(yInterval)));
    } else if(law == 3) {
      core.push_back(Law3(std::move(xInterval), std::move(yInterval)));
    } else if(law == 4) {
      core.push_back(Law4(std::move(xInterval), std::move(yInterval)));
    } else if(law == 5) {
      core.push_back(Law5(std::move(xInterval), std::move(yInterval)));
    } else {
      // TODO unknown interpolation law
      throw std::runtime_error("Unkown interpolation law " + std::to_string(law));
    }

    low = hi - 1;

    // Check for discontinuity
    if(low < ranges::size(x) - 1) {
      if(*(x.begin()+low) == *(x.begin()+low+1)) {low++;}
    }
  }

  Tab1 table(std::move(core));
  
  return table;
}

#endif
