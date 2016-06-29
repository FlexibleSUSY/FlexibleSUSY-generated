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

// File generated at Wed 29 Jun 2016 11:27:09

#include "SM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd SM_input_parameters::get() const
{
   Eigen::ArrayXd pars(3);

   pars(0) = LambdaIN;
   pars(1) = Qin;
   pars(2) = QEWSB;

   return pars;
}

void SM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   LambdaIN = pars(0);
   Qin = pars(1);
   QEWSB = pars(2);

}

std::ostream& operator<<(std::ostream& ostr, const SM_input_parameters& input)
{
   ostr << "LambdaIN = " << INPUT(LambdaIN) << ", ";
   ostr << "Qin = " << INPUT(Qin) << ", ";
   ostr << "QEWSB = " << INPUT(QEWSB) << ", ";

   return ostr;
}

} // namespace flexiblesusy
