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
 * @file MSSM_br_l_to_l_gamma.hpp
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#ifndef MSSM_BrLToLGamma
#define MSSM_BrLToLGamma

namespace softsusy {
class QedQcd;
}

namespace flexiblesusy {
class Physical_input;
class MSSM_mass_eigenstates;

namespace MSSM_br_l_to_l_gamma {
template <typename FIn, typename FOut, typename T1, typename T2>
double lepton_total_decay_width(
      T1 const&, T2 const&,
      const MSSM_mass_eigenstates&, const softsusy::QedQcd&);
double calculate_Fe_to_Fe_VP(
int generationIndex1, int generationIndex2, 
const MSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd, const Physical_input& physical_input);
}
} // namespace flexiblesusy

#endif
