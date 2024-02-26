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
 * @file TMSSM_edm.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef TMSSM_EDM_H
#define TMSSM_EDM_H

#include "TMSSM_mass_eigenstates.hpp"
#include "lowe.h"

namespace flexiblesusy {

class TMSSM_mass_eigenstates;

namespace TMSSM_edm {

template <typename Lepton>
double calculate_edm(const TMSSM_mass_eigenstates&, const softsusy::QedQcd&, int idx);

} // namespace TMSSM_edm
} // namespace flexiblesusy

#endif
