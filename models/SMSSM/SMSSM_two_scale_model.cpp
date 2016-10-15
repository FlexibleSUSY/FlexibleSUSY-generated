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

// File generated at Sat 15 Oct 2016 15:36:24

/**
 * @file SMSSM_two_scale_model.cpp
 * @brief implementation of the SMSSM model class
 *
 * Contains the definition of the SMSSM model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Sat 15 Oct 2016 15:36:24 with FlexibleSUSY
 * 1.7.1 (git commit: 1c1e3234ccd2a3935de013cbdabfb338bedc9204) and SARAH 4.9.1 .
 */

#include "SMSSM_two_scale_model.hpp"

namespace flexiblesusy {

using namespace SMSSM_info;

#define CLASSNAME SMSSM<Two_scale>

CLASSNAME::SMSSM(const SMSSM_input_parameters& input_)
   : Two_scale_model()
   , SMSSM_mass_eigenstates(input_)
{
}

CLASSNAME::~SMSSM()
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
