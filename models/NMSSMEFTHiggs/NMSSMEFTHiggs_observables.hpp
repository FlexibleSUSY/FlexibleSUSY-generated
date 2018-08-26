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

// File generated at Sun 26 Aug 2018 14:04:52

#ifndef NMSSMEFTHiggs_OBSERVABLES_H
#define NMSSMEFTHiggs_OBSERVABLES_H

#include <string>
#include <vector>
#include <Eigen/Core>

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {

class NMSSMEFTHiggs_mass_eigenstates;
class Physical_input;

struct NMSSMEFTHiggs_observables {
   static const int NUMBER_OF_OBSERVABLES = 1;

   NMSSMEFTHiggs_observables();
   Eigen::ArrayXd get() const; ///< returns vector of all observables
   static std::vector<std::string> get_names(); ///< returns vector of all observable names
   void clear(); ///< sets all observables to zero
   void set(const Eigen::ArrayXd&); ///< sets all observables from given vector

   double a_muon; ///< a_muon = (g-2)/2 of the muon (calculated with FlexibleSUSY)

};

NMSSMEFTHiggs_observables calculate_observables(
   NMSSMEFTHiggs_mass_eigenstates&, const softsusy::QedQcd&,
   const Physical_input&);

NMSSMEFTHiggs_observables calculate_observables(
   NMSSMEFTHiggs_mass_eigenstates&, const softsusy::QedQcd&,
   const Physical_input&, double scale);

} // namespace flexiblesusy

#endif
