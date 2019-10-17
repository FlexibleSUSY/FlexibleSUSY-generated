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

// File generated at Wed 16 Oct 2019 22:10:46

/**
 * @file lowNMSSMTanBetaAtMZ_two_scale_model.cpp
 * @brief implementation of the lowNMSSMTanBetaAtMZ model class
 *
 * Contains the definition of the lowNMSSMTanBetaAtMZ model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Wed 16 Oct 2019 22:10:46 with FlexibleSUSY
 * 2.4.1 (git commit: 3e3a10f4fde301d99a4732f29d14f4ac1c646945) and SARAH 4.14.3 .
 */

#include "lowNMSSMTanBetaAtMZ_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME lowNMSSMTanBetaAtMZ<Two_scale>

CLASSNAME::lowNMSSMTanBetaAtMZ(const lowNMSSMTanBetaAtMZ_input_parameters& input_)
   : lowNMSSMTanBetaAtMZ_mass_eigenstates(input_)
{
}

void CLASSNAME::calculate_spectrum()
{
   lowNMSSMTanBetaAtMZ_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   lowNMSSMTanBetaAtMZ_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return lowNMSSMTanBetaAtMZ_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   lowNMSSMTanBetaAtMZ_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   lowNMSSMTanBetaAtMZ_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   lowNMSSMTanBetaAtMZ_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const lowNMSSMTanBetaAtMZ<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
