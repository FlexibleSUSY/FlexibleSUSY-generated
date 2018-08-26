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

// File generated at Sun 26 Aug 2018 15:17:26

/**
 * @file MSSMNoFVatMGUT_a_muon.hpp
 *
 * This file was generated at Sun 26 Aug 2018 15:17:26 with FlexibleSUSY
 * 2.2.0 and SARAH 4.13.0 .
 */

#ifndef MSSMNoFVatMGUT_A_MUON_H
#define MSSMNoFVatMGUT_A_MUON_H

namespace flexiblesusy {
class MSSMNoFVatMGUT_mass_eigenstates;

namespace MSSMNoFVatMGUT_a_muon {
/**
* @fn calculate_a_muon
* @brief Calculates \f$a_\mu = (g-2)_\mu/2\f$ of the muon.
*/
double calculate_a_muon(const MSSMNoFVatMGUT_mass_eigenstates& model);

/**
* @fn calculate_a_muon_uncertainty
* @brief Calculates \f$\Delta a_\mu\f$ of the muon.
*/
double calculate_a_muon_uncertainty(const MSSMNoFVatMGUT_mass_eigenstates& model);
} // namespace MSSMNoFVatMGUT_a_muon
} // namespace flexiblesusy

#endif
