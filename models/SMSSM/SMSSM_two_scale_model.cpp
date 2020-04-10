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

// File generated at Fri 10 Apr 2020 20:19:23

/**
 * @file SMSSM_two_scale_model.cpp
 * @brief implementation of the SMSSM model class
 *
 * Contains the definition of the SMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Fri 10 Apr 2020 20:19:23 with FlexibleSUSY
 * 2.4.2 (git commit: a94199e5620b8684f5d30d0eece5757a5a72c4a4) and SARAH 4.14.3 .
 */

#include "SMSSM_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME SMSSM<Two_scale>

CLASSNAME::SMSSM(const SMSSM_input_parameters& input_)
   : SMSSM_mass_eigenstates(input_)
{
}

void CLASSNAME::calculate_spectrum()
{
   SMSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   SMSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return SMSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   SMSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   SMSSM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   SMSSM_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const SMSSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
