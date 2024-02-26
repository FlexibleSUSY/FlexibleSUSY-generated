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
 * @file THDMII_two_scale_model.cpp
 * @brief implementation of the THDMII model class
 *
 * Contains the definition of the THDMII model class methods
 * which solve EWSB and calculate pole masses and mixings from MSbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.8.0 and SARAH 4.15.1 .
 */

#include "THDMII_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME THDMII<Two_scale>

CLASSNAME::THDMII(const THDMII_slha& model_, bool do_convert_masses_to_slha)
   : THDMII_slha(model_, do_convert_masses_to_slha)
{
}

CLASSNAME::THDMII(const THDMII_input_parameters& input_, bool do_convert_masses_to_slha)
   : THDMII_slha(input_, do_convert_masses_to_slha)
{
}

void CLASSNAME::calculate_spectrum()
{
   THDMII_slha::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   THDMII_slha::clear_problems();
}

std::string CLASSNAME::name() const
{
   return THDMII_slha::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   THDMII_slha::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   THDMII_slha::print(out);
}

void CLASSNAME::set_precision(double p)
{
   THDMII_slha::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const THDMII<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
