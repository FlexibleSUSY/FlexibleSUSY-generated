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
 * @file SM_l_to_l_conversion.hpp
 *
 * This file was generated at Sun 25 Feb 2024 20:21:00 with FlexibleSUSY
 * 2.8.0 and SARAH 4.15.1 .
 */

#ifndef SM_LToLConversion_H
#define SM_LToLConversion_H
#include "lowe.h"

namespace flexiblesusy {
class LToLConversion_settings;
namespace SM_l_to_l_conversion {

enum class Nucleus { Au, Al, Ti };

Eigen::Array<std::complex<double>,13,1> calculate_FeFe_forAll_1loop(int in, int out, const SM_l_to_l_conversion::Nucleus n, const SM_mass_eigenstates& model, const LToLConversion_settings& parameters, const softsusy::QedQcd& qedqcd);

} // namespace SM_l_to_l_conversion
} // namespace flexiblesusy

#endif
