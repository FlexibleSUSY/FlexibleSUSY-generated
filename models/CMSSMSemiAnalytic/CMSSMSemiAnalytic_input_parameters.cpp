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

// File generated at Mon 5 Mar 2018 19:08:55

#include "CMSSMSemiAnalytic_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd CMSSMSemiAnalytic_input_parameters::get() const
{
   Eigen::ArrayXd pars(4);

   pars(0) = m12;
   pars(1) = TanBeta;
   pars(2) = Azero;
   pars(3) = MuInput;

   return pars;
}

void CMSSMSemiAnalytic_input_parameters::set(const Eigen::ArrayXd& pars)
{
   m12 = pars(0);
   TanBeta = pars(1);
   Azero = pars(2);
   MuInput = pars(3);

}

std::ostream& operator<<(std::ostream& ostr, const CMSSMSemiAnalytic_input_parameters& input)
{
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "MuInput = " << INPUT(MuInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
