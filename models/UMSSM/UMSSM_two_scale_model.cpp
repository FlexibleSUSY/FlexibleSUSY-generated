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

// File generated at Sun 10 Jan 2016 15:40:57

/**
 * @file UMSSM_two_scale_model.cpp
 * @brief implementation of the UMSSM model class
 *
 * Contains the definition of the UMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Sun 10 Jan 2016 15:40:57 with FlexibleSUSY
 * 1.3.2 (git commit: v1.3.2) and SARAH 4.6.0 .
 */

#include "UMSSM_two_scale_model.hpp"

namespace flexiblesusy {

using namespace UMSSM_info;

#define CLASSNAME UMSSM<Two_scale>

CLASSNAME::UMSSM(const UMSSM_input_parameters& input_)
   : Two_scale_model()
   , UMSSM_mass_eigenstates(input_)
{
}

CLASSNAME::~UMSSM()
{
}

void CLASSNAME::calculate_spectrum()
{
   UMSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   UMSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return UMSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   UMSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   UMSSM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   UMSSM_mass_eigenstates::set_precision(p);
}

} // namespace flexiblesusy
