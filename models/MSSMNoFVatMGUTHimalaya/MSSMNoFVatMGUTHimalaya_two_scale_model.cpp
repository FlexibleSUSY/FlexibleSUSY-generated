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

// File generated at Mon 5 Mar 2018 14:59:34

/**
 * @file MSSMNoFVatMGUTHimalaya_two_scale_model.cpp
 * @brief implementation of the MSSMNoFVatMGUTHimalaya model class
 *
 * Contains the definition of the MSSMNoFVatMGUTHimalaya model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Mon 5 Mar 2018 14:59:34 with FlexibleSUSY
 * 2.1.0 (git commit: 8f20f6c9c42c159c1588fbc0bb3e15ce5ab6ace3) and SARAH 4.12.3 .
 */

#include "MSSMNoFVatMGUTHimalaya_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME MSSMNoFVatMGUTHimalaya<Two_scale>

CLASSNAME::MSSMNoFVatMGUTHimalaya(const MSSMNoFVatMGUTHimalaya_input_parameters& input_)
   : MSSMNoFVatMGUTHimalaya_mass_eigenstates(input_)
{
}

void CLASSNAME::calculate_spectrum()
{
   MSSMNoFVatMGUTHimalaya_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   MSSMNoFVatMGUTHimalaya_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return MSSMNoFVatMGUTHimalaya_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   MSSMNoFVatMGUTHimalaya_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   MSSMNoFVatMGUTHimalaya_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   MSSMNoFVatMGUTHimalaya_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const MSSMNoFVatMGUTHimalaya<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy