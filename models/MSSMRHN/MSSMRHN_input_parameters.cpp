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

// File generated at Thu 15 Dec 2016 13:00:17

#include "MSSMRHN_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd MSSMRHN_input_parameters::get() const
{
   Eigen::ArrayXd pars(14);

   pars(0) = m0;
   pars(1) = m12;
   pars(2) = TanBeta;
   pars(3) = SignMu;
   pars(4) = Azero;
   pars(5) = BMvInput(0,0);
   pars(6) = BMvInput(0,1);
   pars(7) = BMvInput(0,2);
   pars(8) = BMvInput(1,0);
   pars(9) = BMvInput(1,1);
   pars(10) = BMvInput(1,2);
   pars(11) = BMvInput(2,0);
   pars(12) = BMvInput(2,1);
   pars(13) = BMvInput(2,2);

   return pars;
}

void MSSMRHN_input_parameters::set(const Eigen::ArrayXd& pars)
{
   m0 = pars(0);
   m12 = pars(1);
   TanBeta = pars(2);
   SignMu = pars(3);
   Azero = pars(4);
   BMvInput(0,0) = pars(5);
   BMvInput(0,1) = pars(6);
   BMvInput(0,2) = pars(7);
   BMvInput(1,0) = pars(8);
   BMvInput(1,1) = pars(9);
   BMvInput(1,2) = pars(10);
   BMvInput(2,0) = pars(11);
   BMvInput(2,1) = pars(12);
   BMvInput(2,2) = pars(13);

}

std::ostream& operator<<(std::ostream& ostr, const MSSMRHN_input_parameters& input)
{
   ostr << "m0 = " << INPUT(m0) << ", ";
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "SignMu = " << INPUT(SignMu) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "BMvInput = " << INPUT(BMvInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
