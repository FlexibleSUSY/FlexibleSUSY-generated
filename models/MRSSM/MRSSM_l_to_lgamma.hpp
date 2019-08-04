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

// File generated at Sun 4 Aug 2019 19:30:46

/**
 * @file MRSSM_edm.hpp
 *
 * This file was generated at Sun 4 Aug 2019 19:30:46 with FlexibleSUSY
 * 2.4.0 and SARAH 4.14.2 .
 */

#ifndef MRSSM_MuToEGamma_H
#define MRSSM_MuToEGamma_H

#include "lowe.h"
#include "physical_input.hpp"

namespace flexiblesusy {
class MRSSM_mass_eigenstates;

namespace MRSSM_l_to_lgamma {
template <typename FIn, typename FOut, typename T1, typename T2>
double lepton_total_decay_width(
      T1 const&, T2 const&, 
      const MRSSM_mass_eigenstates&, const softsusy::QedQcd&);

}
} // namespace flexiblesusy

#endif
