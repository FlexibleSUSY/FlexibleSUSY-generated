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
 * @file lowNMSSM_unitarity.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef lowNMSSM_UNITARITY_H
#define lowNMSSM_UNITARITY_H

#include "unitarity.hpp"

namespace flexiblesusy {
class lowNMSSM_mass_eigenstates;

namespace lowNMSSM_unitarity {

UnitarityInfiniteS
max_scattering_eigenvalue_infinite_s(lowNMSSM_mass_eigenstates const&);

} // namespace lowNMSSM_unitarity
} // namespace flexiblesusy

#endif
