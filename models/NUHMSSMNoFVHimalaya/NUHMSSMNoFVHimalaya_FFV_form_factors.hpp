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
 * @file NUHMSSMNoFVHimalaya_FFV_form_factors.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef NUHMSSMNoFVHimalaya_FFVFormFactors_H
#define NUHMSSMNoFVHimalaya_FFVFormFactors_H

#include <valarray>

#include "cxx_qft/NUHMSSMNoFVHimalaya_qft.hpp"

namespace flexiblesusy {

class NUHMSSMNoFVHimalaya_mass_eigenstates;

namespace NUHMSSMNoFVHimalaya_FFV_form_factors {

std::valarray<std::complex<double>> calculate_Fm_Fm_VP_form_factors (
     const NUHMSSMNoFVHimalaya_mass_eigenstates& model, bool discard_SM_contributions);

template <typename Fj, typename Fi, typename V>
std::enable_if_t<
   std::is_same_v<Fj, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm> && std::is_same_v<Fi, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::Fm> && std::is_same_v<V, NUHMSSMNoFVHimalaya_cxx_diagrams::fields::VP>,
   std::valarray<std::complex<double>>
>
calculate_form_factors(     const NUHMSSMNoFVHimalaya_mass_eigenstates&, bool);

} // namespace NUHMSSMNoFVHimalaya_FFV_form_factors
} // namespace flexiblesusy

#endif

