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
 * @file MRSSM2_unitarity.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef MRSSM2_UNITARITY_H
#define MRSSM2_UNITARITY_H

#include "unitarity.hpp"

namespace flexiblesusy {
class MRSSM2_mass_eigenstates;

namespace MRSSM2_unitarity {

UnitarityInfiniteS
max_scattering_eigenvalue_infinite_s(MRSSM2_mass_eigenstates const&);

} // namespace MRSSM2_unitarity
} // namespace flexiblesusy

#endif

