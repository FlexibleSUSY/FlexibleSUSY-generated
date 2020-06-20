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


#include "THDMIIMSSMBC_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd THDMIIMSSMBC_input_parameters::get() const
{
   Eigen::ArrayXd pars(9);

   pars(0) = TanBeta;
   pars(1) = MSUSY;
   pars(2) = MEWSB;
   pars(3) = MuInput;
   pars(4) = MAInput;
   pars(5) = AtInput;
   pars(6) = AbInput;
   pars(7) = AtauInput;
   pars(8) = LambdaLoopOrder;

   return pars;
}

void THDMIIMSSMBC_input_parameters::set(const Eigen::ArrayXd& pars)
{
   TanBeta = pars(0);
   MSUSY = pars(1);
   MEWSB = pars(2);
   MuInput = pars(3);
   MAInput = pars(4);
   AtInput = pars(5);
   AbInput = pars(6);
   AtauInput = pars(7);
   LambdaLoopOrder = pars(8);

}

std::ostream& operator<<(std::ostream& ostr, const THDMIIMSSMBC_input_parameters& input)
{
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "MSUSY = " << INPUT(MSUSY) << ", ";
   ostr << "MEWSB = " << INPUT(MEWSB) << ", ";
   ostr << "MuInput = " << INPUT(MuInput) << ", ";
   ostr << "MAInput = " << INPUT(MAInput) << ", ";
   ostr << "AtInput = " << INPUT(AtInput) << ", ";
   ostr << "AbInput = " << INPUT(AbInput) << ", ";
   ostr << "AtauInput = " << INPUT(AtauInput) << ", ";
   ostr << "LambdaLoopOrder = " << INPUT(LambdaLoopOrder) << ", ";

   return ostr;
}

} // namespace flexiblesusy
