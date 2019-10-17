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
 * @file CE6SSM_FFV_form_factors.hpp
 *
 * This file was generated at Wed 16 Oct 2019 20:12:36 with FlexibleSUSY
 * 2.4.1 and SARAH 4.14.3 .
 */

#ifndef CE6SSM_FFVFormFactors_H
#define CE6SSM_FFVFormFactors_H

#include <valarray>

#include "cxx_qft/CE6SSM_qft.hpp"

namespace flexiblesusy {

class CE6SSM_mass_eigenstates;

namespace CE6SSM_FFV_form_factors {
std::valarray<std::complex<double>> calculate_Fe_Fe_VP_form_factors (
   int generationIndex1, int generationIndex2,
   const CE6SSM_mass_eigenstates& model, bool discard_SM_contributions);
}

} // namespace flexiblesusy

#endif
