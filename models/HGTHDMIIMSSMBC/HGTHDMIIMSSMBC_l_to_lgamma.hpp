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

// File generated at Fri 10 Apr 2020 19:28:43

/**
 * @file HGTHDMIIMSSMBC_edm.hpp
 *
 * This file was generated at Fri 10 Apr 2020 19:28:43 with FlexibleSUSY
 * 2.4.2 and SARAH 4.14.3 .
 */

#ifndef HGTHDMIIMSSMBC_MuToEGamma_H
#define HGTHDMIIMSSMBC_MuToEGamma_H

#include "lowe.h"
#include "physical_input.hpp"

namespace flexiblesusy {
class HGTHDMIIMSSMBC_mass_eigenstates;

namespace HGTHDMIIMSSMBC_l_to_lgamma {
template <typename FIn, typename FOut, typename T1, typename T2>
double lepton_total_decay_width(
      T1 const&, T2 const&, 
      const HGTHDMIIMSSMBC_mass_eigenstates&, const softsusy::QedQcd&);

}
} // namespace flexiblesusy

#endif
