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

// File generated at Tue 8 Mar 2016 18:44:38

/**
 * @file MSSM_two_scale_model.cpp
 * @brief implementation of the MSSM model class
 *
 * Contains the definition of the MSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Tue 8 Mar 2016 18:44:38 with FlexibleSUSY
 * 1.4.0 (git commit: v1.4.0) and SARAH 4.7.0 .
 */

#include "MSSM_two_scale_model.hpp"

namespace flexiblesusy {

using namespace MSSM_info;

#define CLASSNAME MSSM<Two_scale>

CLASSNAME::MSSM(const MSSM_input_parameters& input_)
   : Two_scale_model()
   , MSSM_mass_eigenstates(input_)
{
}

CLASSNAME::~MSSM()
{
}

void CLASSNAME::calculate_spectrum()
{
   MSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   MSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return MSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   MSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   MSSM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   MSSM_mass_eigenstates::set_precision(p);
}

} // namespace flexiblesusy
