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

// File generated at Wed 16 Oct 2019 19:59:26

/**
 * @file MSSMEFTHiggs_a_muon.hpp
 *
 * This file was generated at Wed 16 Oct 2019 19:59:26 with FlexibleSUSY
 * 2.4.1 and SARAH 4.14.3 .
 */

#ifndef MSSMEFTHiggs_A_MUON_H
#define MSSMEFTHiggs_A_MUON_H

namespace softsusy {
   class QedQcd;
} // namespace softsusy

namespace flexiblesusy {
class MSSMEFTHiggs_mass_eigenstates;

namespace MSSMEFTHiggs_a_muon {
/**
* @fn calculate_a_muon
* @brief Calculates \f$a_\mu = (g-2)_\mu/2\f$ of the muon.
*/
double calculate_a_muon(const MSSMEFTHiggs_mass_eigenstates& model, const softsusy::QedQcd& qedqcd);

/**
* @fn calculate_a_muon_uncertainty
* @brief Calculates \f$\Delta a_\mu\f$ of the muon.
*/
double calculate_a_muon_uncertainty(const MSSMEFTHiggs_mass_eigenstates& model, const softsusy::QedQcd& qedqcd);
} // namespace MSSMEFTHiggs_a_muon
} // namespace flexiblesusy

#endif
