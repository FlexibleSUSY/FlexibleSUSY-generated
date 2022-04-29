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
 * @file MSSMNoFVHimalaya_two_scale_model.cpp
 * @brief implementation of the MSSMNoFVHimalaya model class
 *
 * Contains the definition of the MSSMNoFVHimalaya model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.6.2 and SARAH 4.14.5 .
 */

#include "MSSMNoFVHimalaya_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME MSSMNoFVHimalaya<Two_scale>

CLASSNAME::MSSMNoFVHimalaya(const MSSMNoFVHimalaya_slha& model_, bool do_convert_masses_to_slha)
   : MSSMNoFVHimalaya_slha(model_, do_convert_masses_to_slha)
{
}

CLASSNAME::MSSMNoFVHimalaya(const MSSMNoFVHimalaya_input_parameters& input_, bool do_convert_masses_to_slha)
   : MSSMNoFVHimalaya_slha(input_, do_convert_masses_to_slha)
{
}

void CLASSNAME::calculate_spectrum()
{
   MSSMNoFVHimalaya_slha::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   MSSMNoFVHimalaya_slha::clear_problems();
}

std::string CLASSNAME::name() const
{
   return MSSMNoFVHimalaya_slha::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   MSSMNoFVHimalaya_slha::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   MSSMNoFVHimalaya_slha::print(out);
}

void CLASSNAME::set_precision(double p)
{
   MSSMNoFVHimalaya_slha::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFVHimalaya<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
