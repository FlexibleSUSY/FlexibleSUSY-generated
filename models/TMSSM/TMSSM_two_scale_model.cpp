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

// File generated at Sun 26 Aug 2018 14:15:53

/**
 * @file TMSSM_two_scale_model.cpp
 * @brief implementation of the TMSSM model class
 *
 * Contains the definition of the TMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Sun 26 Aug 2018 14:15:53 with FlexibleSUSY
 * 2.2.0 (git commit: 8489097de2d6938a6da0149378457b5ad13d9425) and SARAH 4.13.0 .
 */

#include "TMSSM_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME TMSSM<Two_scale>

CLASSNAME::TMSSM(const TMSSM_input_parameters& input_)
   : TMSSM_mass_eigenstates(input_)
{
}

void CLASSNAME::calculate_spectrum()
{
   TMSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   TMSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return TMSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   TMSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   TMSSM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   TMSSM_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const TMSSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
