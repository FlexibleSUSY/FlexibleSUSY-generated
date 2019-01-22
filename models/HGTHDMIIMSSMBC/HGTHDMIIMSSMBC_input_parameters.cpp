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

// File generated at Tue 22 Jan 2019 16:22:38

#include "HGTHDMIIMSSMBC_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd HGTHDMIIMSSMBC_input_parameters::get() const
{
   Eigen::ArrayXd pars(12);

   pars(0) = TanBeta;
   pars(1) = MSUSY;
   pars(2) = MEWSB;
   pars(3) = MuInput;
   pars(4) = M1Input;
   pars(5) = M2Input;
   pars(6) = M3Input;
   pars(7) = MAInput;
   pars(8) = AtInput;
   pars(9) = AbInput;
   pars(10) = AtauInput;
   pars(11) = LambdaLoopOrder;

   return pars;
}

void HGTHDMIIMSSMBC_input_parameters::set(const Eigen::ArrayXd& pars)
{
   TanBeta = pars(0);
   MSUSY = pars(1);
   MEWSB = pars(2);
   MuInput = pars(3);
   M1Input = pars(4);
   M2Input = pars(5);
   M3Input = pars(6);
   MAInput = pars(7);
   AtInput = pars(8);
   AbInput = pars(9);
   AtauInput = pars(10);
   LambdaLoopOrder = pars(11);

}

std::ostream& operator<<(std::ostream& ostr, const HGTHDMIIMSSMBC_input_parameters& input)
{
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "MSUSY = " << INPUT(MSUSY) << ", ";
   ostr << "MEWSB = " << INPUT(MEWSB) << ", ";
   ostr << "MuInput = " << INPUT(MuInput) << ", ";
   ostr << "M1Input = " << INPUT(M1Input) << ", ";
   ostr << "M2Input = " << INPUT(M2Input) << ", ";
   ostr << "M3Input = " << INPUT(M3Input) << ", ";
   ostr << "MAInput = " << INPUT(MAInput) << ", ";
   ostr << "AtInput = " << INPUT(AtInput) << ", ";
   ostr << "AbInput = " << INPUT(AbInput) << ", ";
   ostr << "AtauInput = " << INPUT(AtauInput) << ", ";
   ostr << "LambdaLoopOrder = " << INPUT(LambdaLoopOrder) << ", ";

   return ostr;
}

} // namespace flexiblesusy
