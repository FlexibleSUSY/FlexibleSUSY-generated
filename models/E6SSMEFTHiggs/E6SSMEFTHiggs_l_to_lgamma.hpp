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

// File generated at Wed 16 Oct 2019 20:19:34

/**
 * @file E6SSMEFTHiggs_edm.hpp
 *
 * This file was generated at Wed 16 Oct 2019 20:19:34 with FlexibleSUSY
 * 2.4.1 and SARAH 4.14.3 .
 */

#ifndef E6SSMEFTHiggs_MuToEGamma_H
#define E6SSMEFTHiggs_MuToEGamma_H

#include "lowe.h"
#include "physical_input.hpp"

namespace flexiblesusy {
class E6SSMEFTHiggs_mass_eigenstates;

namespace E6SSMEFTHiggs_l_to_lgamma {
template <typename FIn, typename FOut, typename T1, typename T2>
double lepton_total_decay_width(
      T1 const&, T2 const&, 
      const E6SSMEFTHiggs_mass_eigenstates&, const softsusy::QedQcd&);

}
} // namespace flexiblesusy

#endif
