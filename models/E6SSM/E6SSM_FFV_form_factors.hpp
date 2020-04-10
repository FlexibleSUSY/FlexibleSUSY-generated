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
 * @file E6SSM_FFV_form_factors.hpp
 *
 * This file was generated at Fri 10 Apr 2020 20:05:52 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#ifndef E6SSM_FFVFormFactors_H
#define E6SSM_FFVFormFactors_H

#include <valarray>

#include "cxx_qft/E6SSM_qft.hpp"

namespace flexiblesusy {

class E6SSM_mass_eigenstates;

namespace E6SSM_FFV_form_factors {
std::valarray<std::complex<double>> calculate_Fe_Fe_VP_form_factors (
   int generationIndex1, int generationIndex2,
   const E6SSM_mass_eigenstates& model, bool discard_SM_contributions);
}

} // namespace flexiblesusy

#endif

