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

// File generated at Mon 9 May 2016 13:26:52

/**
 * @file MSSMRHN_two_scale_model.cpp
 * @brief implementation of the MSSMRHN model class
 *
 * Contains the definition of the MSSMRHN model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Mon 9 May 2016 13:26:52 with FlexibleSUSY
 * 1.4.2 (git commit: ba53b7080ae303fc6b5ef4b4ce12d05fef5b6211) and SARAH 4.8.5 .
 */

#include "MSSMRHN_two_scale_model.hpp"

namespace flexiblesusy {

using namespace MSSMRHN_info;

#define CLASSNAME MSSMRHN<Two_scale>

CLASSNAME::MSSMRHN(const MSSMRHN_input_parameters& input_)
   : Two_scale_model()
   , MSSMRHN_mass_eigenstates(input_)
{
}

CLASSNAME::~MSSMRHN()
{
}

void CLASSNAME::calculate_spectrum()
{
   MSSMRHN_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   MSSMRHN_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return MSSMRHN_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   MSSMRHN_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   MSSMRHN_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   MSSMRHN_mass_eigenstates::set_precision(p);
}

} // namespace flexiblesusy
