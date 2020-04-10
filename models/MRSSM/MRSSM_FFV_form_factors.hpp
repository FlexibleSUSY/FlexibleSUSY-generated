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
 * @file MRSSM_FFV_form_factors.hpp
 *
 * This file was generated at Fri 10 Apr 2020 19:59:01 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#ifndef MRSSM_FFVFormFactors_H
#define MRSSM_FFVFormFactors_H

#include <valarray>

#include "cxx_qft/MRSSM_qft.hpp"

namespace flexiblesusy {

class MRSSM_mass_eigenstates;

namespace MRSSM_FFV_form_factors {
std::valarray<std::complex<double>> calculate_Fe_Fe_VP_form_factors (
   int generationIndex1, int generationIndex2,
   const MRSSM_mass_eigenstates& model, bool discard_SM_contributions);
}

} // namespace flexiblesusy

#endif

