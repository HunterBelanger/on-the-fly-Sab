#ifndef TAB1_H
#define TAB1_H

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

#endif