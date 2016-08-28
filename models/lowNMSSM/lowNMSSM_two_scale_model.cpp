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

// File generated at Sun 28 Aug 2016 15:14:30

/**
 * @file lowNMSSM_two_scale_model.cpp
 * @brief implementation of the lowNMSSM model class
 *
 * Contains the definition of the lowNMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Sun 28 Aug 2016 15:14:30 with FlexibleSUSY
 * 1.6.1 (git commit: de78c399bea838d002fe54cbc91eb6f88fa62626) and SARAH 4.9.1 .
 */

#include "lowNMSSM_two_scale_model.hpp"

namespace flexiblesusy {

using namespace lowNMSSM_info;

#define CLASSNAME lowNMSSM<Two_scale>

CLASSNAME::lowNMSSM(const lowNMSSM_input_parameters& input_)
   : Two_scale_model()
   , lowNMSSM_mass_eigenstates(input_)
{
}

CLASSNAME::~lowNMSSM()
{
}

void CLASSNAME::calculate_spectrum()
{
   lowNMSSM_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   lowNMSSM_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return lowNMSSM_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   lowNMSSM_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   lowNMSSM_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   lowNMSSM_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const lowNMSSM<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
