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

// File generated at Thu 15 Dec 2016 12:58:40

/**
 * @file NUTSMSSM_two_scale_model.cpp
 * @brief implementation of the NUTSMSSM model class
 *
 * Contains the definition of the NUTSMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Thu 15 Dec 2016 12:58:40 with FlexibleSUSY
 * 1.7.2 (git commit: 0d19299fef514160cb7541a03abb9b2c3365f927) and SARAH 4.9.1 .
 */

#include "NUTSMSSM_two_scale_model.hpp"

namespace flexiblesusy {

using namespace NUTSMSSM_info;

#define CLASSNAME NUTSMSSM<Two_scale>

CLASSNAME::NUTSMSSM(const NUTSMSSM_input_parameters& input_)
   : Two_scale_model()
   , NUTSMSSM_mass_eigenstates(input_)
{
}

CLASSNAME::~NUTSMSSM()
{
}

void CLASSNAME::calculate_spectrum()
{
   NUTSMSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   NUTSMSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return NUTSMSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   NUTSMSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   NUTSMSSM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   NUTSMSSM_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const NUTSMSSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
