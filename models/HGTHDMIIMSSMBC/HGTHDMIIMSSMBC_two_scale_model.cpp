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

// File generated at Mon 9 May 2016 12:07:06

/**
 * @file HGTHDMIIMSSMBC_two_scale_model.cpp
 * @brief implementation of the HGTHDMIIMSSMBC model class
 *
 * Contains the definition of the HGTHDMIIMSSMBC model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Mon 9 May 2016 12:07:06 with FlexibleSUSY
 * 1.4.2 (git commit: ba53b7080ae303fc6b5ef4b4ce12d05fef5b6211) and SARAH 4.8.5 .
 */

#include "HGTHDMIIMSSMBC_two_scale_model.hpp"

namespace flexiblesusy {

using namespace HGTHDMIIMSSMBC_info;

#define CLASSNAME HGTHDMIIMSSMBC<Two_scale>

CLASSNAME::HGTHDMIIMSSMBC(const HGTHDMIIMSSMBC_input_parameters& input_)
   : Two_scale_model()
   , HGTHDMIIMSSMBC_mass_eigenstates(input_)
{
}

CLASSNAME::~HGTHDMIIMSSMBC()
{
}

void CLASSNAME::calculate_spectrum()
{
   HGTHDMIIMSSMBC_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   HGTHDMIIMSSMBC_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return HGTHDMIIMSSMBC_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   HGTHDMIIMSSMBC_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   HGTHDMIIMSSMBC_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   HGTHDMIIMSSMBC_mass_eigenstates::set_precision(p);
}

} // namespace flexiblesusy