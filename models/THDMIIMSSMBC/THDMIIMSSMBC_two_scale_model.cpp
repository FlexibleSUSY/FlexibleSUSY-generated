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

// File generated at Mon 27 Feb 2017 13:24:46

/**
 * @file THDMIIMSSMBC_two_scale_model.cpp
 * @brief implementation of the THDMIIMSSMBC model class
 *
 * Contains the definition of the THDMIIMSSMBC model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Mon 27 Feb 2017 13:24:46 with FlexibleSUSY
 * 1.7.3 (git commit: 622a80d5da461a0a259a094325cd734ff8e79c61) and SARAH 4.9.3 .
 */

#include "THDMIIMSSMBC_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME THDMIIMSSMBC<Two_scale>

CLASSNAME::THDMIIMSSMBC(const THDMIIMSSMBC_input_parameters& input_)
   : Two_scale_model()
   , THDMIIMSSMBC_mass_eigenstates(input_)
{
}

CLASSNAME::~THDMIIMSSMBC()
{
}

void CLASSNAME::calculate_spectrum()
{
   THDMIIMSSMBC_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   THDMIIMSSMBC_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return THDMIIMSSMBC_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   THDMIIMSSMBC_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   THDMIIMSSMBC_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   THDMIIMSSMBC_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const THDMIIMSSMBC<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
