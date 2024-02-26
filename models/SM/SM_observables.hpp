// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================


#ifndef SM_OBSERVABLES_H
#define SM_OBSERVABLES_H

#include "observable_problems.hpp"
#include "spectrum_generator_settings.hpp"
#include <string>
#include <vector>
#include <Eigen/Core>

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {

class SM_mass_eigenstates;
class Physical_input;
class LToLConversion_settings;

struct SM_observables {
   static constexpr int NUMBER_OF_OBSERVABLES = 55;

   SM_observables();
   Eigen::ArrayXd get() const; ///< returns vector of all observables
   static std::vector<std::string> get_names(); ///< returns vector of all observable names
   void clear(); ///< sets all observables to zero
   void set(const Eigen::ArrayXd&); ///< sets all observables from given vector

   Observable_problems problems;
   double amm_Fe_0; ///< Delta(g-2)/2 of Fe(1) (calculated with FlexibleSUSY)
   double amm_Fe_1; ///< Delta(g-2)/2 of Fe(2) (calculated with FlexibleSUSY)
   double amm_Fe_2; ///< Delta(g-2)/2 of Fe(3) (calculated with FlexibleSUSY)
   Eigen::Array<std::complex<double>,13,1> Fe2to11bar1_All_1loop; ///< BR(Fe(2)->Fe(1)Fe(1)barFe(1)) for All at 1 loop
   Eigen::Array<std::complex<double>,13,1> Fe2Fe1inAl_All_1loop; ///< CR(Fe(2)->Fe(1)) in Al for All at 1 loop

};

SM_observables calculate_observables(
   const SM_mass_eigenstates&,
   const softsusy::QedQcd&,
   const LToLConversion_settings& ltolconversion_settings,
   const Physical_input&,
   const Spectrum_generator_settings&);

SM_observables calculate_observables(
   const SM_mass_eigenstates&,
   const softsusy::QedQcd&,
   const LToLConversion_settings& ltolconversion_settings,
   const Physical_input&,
   const Spectrum_generator_settings&,
   double scale);

} // namespace flexiblesusy

#endif
