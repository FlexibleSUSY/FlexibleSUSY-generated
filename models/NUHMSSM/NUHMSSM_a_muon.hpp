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


/**
 * @file NUHMSSM_a_muon.hpp
 *
 * This file was generated with FlexibleSUSY 2.6.1 and SARAH 4.14.5 .
 */

#ifndef NUHMSSM_A_MUON_H
#define NUHMSSM_A_MUON_H

namespace softsusy {
   class QedQcd;
} // namespace softsusy

namespace flexiblesusy {
class NUHMSSM_mass_eigenstates;

namespace NUHMSSM_a_muon {
/**
* @fn calculate_a_muon
* @brief Calculates \f$a_\mu = (g-2)_\mu/2\f$ of the muon.
*/
double calculate_a_muon(const NUHMSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd);

/**
* @fn calculate_a_muon_uncertainty
* @brief Calculates \f$\Delta a_\mu\f$ of the muon.
*/
double calculate_a_muon_uncertainty(const NUHMSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd);
} // namespace NUHMSSM_a_muon
} // namespace flexiblesusy

#endif
