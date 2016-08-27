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

// File generated at Sat 27 Aug 2016 12:54:57

/**
 * @file NUTNMSSM_two_scale_model.cpp
 * @brief implementation of the NUTNMSSM model class
 *
 * Contains the definition of the NUTNMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Sat 27 Aug 2016 12:54:57 with FlexibleSUSY
 * 1.6.0 (git commit: b48da3168d7e9ce639e93a4ea0b216e3468d6dc3) and SARAH 4.9.1 .
 */

#include "NUTNMSSM_two_scale_model.hpp"

namespace flexiblesusy {

using namespace NUTNMSSM_info;

#define CLASSNAME NUTNMSSM<Two_scale>

CLASSNAME::NUTNMSSM(const NUTNMSSM_input_parameters& input_)
   : Two_scale_model()
   , NUTNMSSM_mass_eigenstates(input_)
{
}

CLASSNAME::~NUTNMSSM()
{
}

void CLASSNAME::calculate_spectrum()
{
   NUTNMSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   NUTNMSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return NUTNMSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   NUTNMSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   NUTNMSSM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   NUTNMSSM_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const NUTNMSSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
