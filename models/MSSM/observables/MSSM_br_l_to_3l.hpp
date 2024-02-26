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
 * @file MSSM_br_l_to_3l.hpp
 *
 * This file was generated at Sun 25 Feb 2024 21:19:20 with FlexibleSUSY
 * 2.8.0 and SARAH 4.15.1 .
 */

#ifndef MSSM_BrLTo3L_H
#define MSSM_BrLTo3L_H
#include "lowe.h"

namespace flexiblesusy {
namespace MSSM_br_l_to_3l {

Eigen::Array<std::complex<double>,13,1> calculate_Fe_to_FeFebarFe_for_All_1loop(int nI, int nO, int nA, const MSSM_mass_eigenstates& model, const softsusy::QedQcd& qedqcd);

} // namespace MSSM_br_l_to_3l
} // namespace flexiblesusy

#endif
