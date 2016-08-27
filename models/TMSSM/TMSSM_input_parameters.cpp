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

// File generated at Sat 27 Aug 2016 11:55:22

#include "TMSSM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd TMSSM_input_parameters::get() const
{
   Eigen::ArrayXd pars(9);

   pars(0) = m0;
   pars(1) = m12;
   pars(2) = TanBeta;
   pars(3) = SignMu;
   pars(4) = Azero;
   pars(5) = MTinput;
   pars(6) = Qin;
   pars(7) = LambdaInput;
   pars(8) = vTInput;

   return pars;
}

void TMSSM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   m0 = pars(0);
   m12 = pars(1);
   TanBeta = pars(2);
   SignMu = pars(3);
   Azero = pars(4);
   MTinput = pars(5);
   Qin = pars(6);
   LambdaInput = pars(7);
   vTInput = pars(8);

}

std::ostream& operator<<(std::ostream& ostr, const TMSSM_input_parameters& input)
{
   ostr << "m0 = " << INPUT(m0) << ", ";
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "SignMu = " << INPUT(SignMu) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "MTinput = " << INPUT(MTinput) << ", ";
   ostr << "Qin = " << INPUT(Qin) << ", ";
   ostr << "LambdaInput = " << INPUT(LambdaInput) << ", ";
   ostr << "vTInput = " << INPUT(vTInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
