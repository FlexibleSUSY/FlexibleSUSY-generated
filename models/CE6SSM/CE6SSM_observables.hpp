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


#ifndef CE6SSM_OBSERVABLES_H
#define CE6SSM_OBSERVABLES_H

#include "observable_problems.hpp"
#include "spectrum_generator_settings.hpp"
#include <string>
#include <vector>
#include <Eigen/Core>

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {

class CE6SSM_mass_eigenstates;
class Physical_input;


struct CE6SSM_observables {
   static constexpr int NUMBER_OF_OBSERVABLES = 1;

   CE6SSM_observables();
   Eigen::ArrayXd get() const; ///< returns vector of all observables
   static std::vector<std::string> get_names(); ///< returns vector of all observable names
   void clear(); ///< sets all observables to zero
   void set(const Eigen::ArrayXd&); ///< sets all observables from given vector

   Observable_problems problems;
   double amm_Fe_1; ///< Delta(g-2)/2 of Fe(2) (calculated with FlexibleSUSY)

};

CE6SSM_observables calculate_observables(
   const CE6SSM_mass_eigenstates&,
   const softsusy::QedQcd&,
   
   const Physical_input&,
   const Spectrum_generator_settings&);

CE6SSM_observables calculate_observables(
   const CE6SSM_mass_eigenstates&,
   const softsusy::QedQcd&,
   
   const Physical_input&,
   const Spectrum_generator_settings&,
   double scale);

} // namespace flexiblesusy

#endif
