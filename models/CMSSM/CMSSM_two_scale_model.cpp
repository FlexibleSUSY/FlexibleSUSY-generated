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

// File generated at Fri 8 Jan 2016 15:32:19

/**
 * @file CMSSM_two_scale_model.cpp
 * @brief implementation of the CMSSM model class
 *
 * Contains the definition of the CMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Fri 8 Jan 2016 15:32:19 with FlexibleSUSY
 * 1.3.1 (git commit: v1.3.1) and SARAH 4.6.0 .
 */

#include "CMSSM_two_scale_model.hpp"

namespace flexiblesusy {

using namespace CMSSM_info;

#define CLASSNAME CMSSM<Two_scale>

CLASSNAME::CMSSM(const CMSSM_input_parameters& input_)
   : Two_scale_model()
   , CMSSM_mass_eigenstates(input_)
{
}

CLASSNAME::~CMSSM()
{
}

void CLASSNAME::calculate_spectrum()
{
   CMSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   CMSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return CMSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   CMSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   CMSSM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   CMSSM_mass_eigenstates::set_precision(p);
}

} // namespace flexiblesusy
