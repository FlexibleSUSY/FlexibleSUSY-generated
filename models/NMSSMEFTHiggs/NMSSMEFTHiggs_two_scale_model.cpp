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
 * @file NMSSMEFTHiggs_two_scale_model.cpp
 * @brief implementation of the NMSSMEFTHiggs model class
 *
 * Contains the definition of the NMSSMEFTHiggs model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated with FlexibleSUSY 2.6.0 and SARAH 4.14.4 .
 */

#include "NMSSMEFTHiggs_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME NMSSMEFTHiggs<Two_scale>

CLASSNAME::NMSSMEFTHiggs(const NMSSMEFTHiggs_slha& model_, bool do_convert_masses_to_slha)
   : NMSSMEFTHiggs_slha(model_, do_convert_masses_to_slha)
{
}

CLASSNAME::NMSSMEFTHiggs(const NMSSMEFTHiggs_input_parameters& input_, bool do_convert_masses_to_slha)
   : NMSSMEFTHiggs_slha(input_, do_convert_masses_to_slha)
{
}

void CLASSNAME::calculate_spectrum()
{
   NMSSMEFTHiggs_slha::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   NMSSMEFTHiggs_slha::clear_problems();
}

std::string CLASSNAME::name() const
{
   return NMSSMEFTHiggs_slha::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   NMSSMEFTHiggs_slha::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   NMSSMEFTHiggs_slha::print(out);
}

void CLASSNAME::set_precision(double p)
{
   NMSSMEFTHiggs_slha::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const NMSSMEFTHiggs<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
