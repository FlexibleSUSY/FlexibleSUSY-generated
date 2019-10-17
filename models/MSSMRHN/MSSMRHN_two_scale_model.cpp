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

// File generated at Wed 16 Oct 2019 22:46:48

/**
 * @file MSSMRHN_two_scale_model.cpp
 * @brief implementation of the MSSMRHN model class
 *
 * Contains the definition of the MSSMRHN model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Wed 16 Oct 2019 22:46:48 with FlexibleSUSY
 * 2.4.1 (git commit: 3e3a10f4fde301d99a4732f29d14f4ac1c646945) and SARAH 4.14.3 .
 */

#include "MSSMRHN_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME MSSMRHN<Two_scale>

CLASSNAME::MSSMRHN(const MSSMRHN_input_parameters& input_)
   : MSSMRHN_mass_eigenstates(input_)
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

std::ostream& operator<<(std::ostream& ostr, const MSSMRHN<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
