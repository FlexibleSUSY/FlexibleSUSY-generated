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
 * @file MRSSMEFTHiggs_FFV_form_factors.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef MRSSMEFTHiggs_FFVFormFactors_H
#define MRSSMEFTHiggs_FFVFormFactors_H

#include <valarray>

#include "cxx_qft/MRSSMEFTHiggs_qft.hpp"

namespace flexiblesusy {

class MRSSMEFTHiggs_mass_eigenstates;

namespace MRSSMEFTHiggs_FFV_form_factors {

std::valarray<std::complex<double>> calculate_Fe_Fe_VP_form_factors (
   int generationIndex1, int generationIndex2,
   const MRSSMEFTHiggs_mass_eigenstates& model, bool discard_SM_contributions);

template <typename Fj, typename Fi, typename V>
std::enable_if_t<
   std::is_same_v<Fj, MRSSMEFTHiggs_cxx_diagrams::fields::Fe> && std::is_same_v<Fi, MRSSMEFTHiggs_cxx_diagrams::fields::Fe> && std::is_same_v<V, MRSSMEFTHiggs_cxx_diagrams::fields::VP>,
   std::valarray<std::complex<double>>
>
calculate_form_factors(   int, int,
const MRSSMEFTHiggs_mass_eigenstates&, bool);

} // namespace MRSSMEFTHiggs_FFV_form_factors
} // namespace flexiblesusy

#endif

