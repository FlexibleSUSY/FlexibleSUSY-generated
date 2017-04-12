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

// File generated at Wed 12 Apr 2017 11:44:42

/**
 * @file MRSSMtower_two_scale_model.cpp
 * @brief implementation of the MRSSMtower model class
 *
 * Contains the definition of the MRSSMtower model class methods
 * which solve EWSB and calculate pole masses and mixings from DRbar
 * parameters.
 *
 * This file was generated at Wed 12 Apr 2017 11:44:42 with FlexibleSUSY
 * 1.7.4 (git commit: bf9e92a2ddb43c203483621f6150a96a16f51536) and SARAH 4.11.0 .
 */

#include "MRSSMtower_two_scale_model.hpp"

namespace flexiblesusy {

#define CLASSNAME MRSSMtower<Two_scale>

CLASSNAME::MRSSMtower(const MRSSMtower_input_parameters& input_)
   : Two_scale_model()
   , MRSSMtower_mass_eigenstates(input_)
{
}

CLASSNAME::~MRSSMtower()
{
}

void CLASSNAME::calculate_spectrum()
{
   MRSSMtower_mass_eigenstates::calculate_spectrum();
}

void CLASSNAME::clear_problems()
{
   MRSSMtower_mass_eigenstates::clear_problems();
}

std::string CLASSNAME::name() const
{
   return MRSSMtower_mass_eigenstates::name();
}

void CLASSNAME::run_to(double scale, double eps)
{
   MRSSMtower_mass_eigenstates::run_to(scale, eps);
}

void CLASSNAME::print(std::ostream& out) const
{
   MRSSMtower_mass_eigenstates::print(out);
}

void CLASSNAME::set_precision(double p)
{
   MRSSMtower_mass_eigenstates::set_precision(p);
}

std::ostream& operator<<(std::ostream& ostr, const MRSSMtower<Two_scale>& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
