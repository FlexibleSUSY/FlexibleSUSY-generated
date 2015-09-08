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

// File generated at Tue 8 Sep 2015 13:57:18

/**
 * @file CMSSMNoFV_two_scale_model.cpp
 * @brief implementation of the CMSSMNoFV model class
 *
 * Contains the definition of the CMSSMNoFV model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Tue 8 Sep 2015 13:57:18 with FlexibleSUSY
 * 1.2.2 (git commit: v1.2.2) and SARAH 4.5.8 .
 */

#include "CMSSMNoFV_two_scale_model.hpp"

namespace flexiblesusy {

using namespace CMSSMNoFV_info;

#define CLASSNAME CMSSMNoFV<Two_scale>

CLASSNAME::CMSSMNoFV(const CMSSMNoFV_input_parameters& input_)
   : Two_scale_model()
   , CMSSMNoFV_mass_eigenstates(input_)
{
}

CLASSNAME::~CMSSMNoFV()
{
}

void CLASSNAME::calculate_spectrum()
{
   CMSSMNoFV_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   CMSSMNoFV_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return CMSSMNoFV_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   CMSSMNoFV_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   CMSSMNoFV_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   CMSSMNoFV_mass_eigenstates::set_precision(p);
}

} // namespace flexiblesusy
