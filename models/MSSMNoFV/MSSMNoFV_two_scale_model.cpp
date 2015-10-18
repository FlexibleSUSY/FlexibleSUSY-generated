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

// File generated at Sun 18 Oct 2015 13:49:15

/**
 * @file MSSMNoFV_two_scale_model.cpp
 * @brief implementation of the MSSMNoFV model class
 *
 * Contains the definition of the MSSMNoFV model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Sun 18 Oct 2015 13:49:15 with FlexibleSUSY
 * 1.2.3 (git commit: v1.2.3-2-g5f1c55e) and SARAH 4.5.8 .
 */

#include "MSSMNoFV_two_scale_model.hpp"

namespace flexiblesusy {

using namespace MSSMNoFV_info;

#define CLASSNAME MSSMNoFV<Two_scale>

CLASSNAME::MSSMNoFV(const MSSMNoFV_input_parameters& input_)
   : Two_scale_model()
   , MSSMNoFV_mass_eigenstates(input_)
{
}

CLASSNAME::~MSSMNoFV()
{
}

void CLASSNAME::calculate_spectrum()
{
   MSSMNoFV_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   MSSMNoFV_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return MSSMNoFV_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   MSSMNoFV_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   MSSMNoFV_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   MSSMNoFV_mass_eigenstates::set_precision(p);
}

} // namespace flexiblesusy
