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
 * @file lowMSSM_amm.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef lowMSSM_AMM_H
#define lowMSSM_AMM_H

namespace softsusy {
   class QedQcd;
} // namespace softsusy

namespace flexiblesusy {
class Spectrum_generator_settings;

class lowMSSM_mass_eigenstates;

namespace lowMSSM_amm {
/**
* @fn calculate_amm
* @brief Calculates \f$a = Delta(g-2)/2\f$ of the lepton.
*/
template <typename Lepton>
double calculate_amm(const lowMSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings&, int idx);

/**
* @fn calculate_amm_uncertainty
* @brief Calculates uncertainty of \f$\Delta a\f$ of the lepton.
*/
template <typename Lepton>
double calculate_amm_uncertainty(const lowMSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Spectrum_generator_settings&, int idx);
} // namespace lowMSSM_amm
} // namespace flexiblesusy

#endif
