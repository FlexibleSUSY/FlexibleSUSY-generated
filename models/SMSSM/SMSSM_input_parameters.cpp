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

// File generated at Sun 26 Aug 2018 14:31:55

#include "SMSSM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd SMSSM_input_parameters::get() const
{
   Eigen::ArrayXd pars(11);

   pars(0) = m0;
   pars(1) = m12;
   pars(2) = TanBeta;
   pars(3) = SignMu;
   pars(4) = Azero;
   pars(5) = LambdaInput;
   pars(6) = KappaInput;
   pars(7) = LambdaSInput;
   pars(8) = L1Input;
   pars(9) = MSInput;
   pars(10) = BMSInput;

   return pars;
}

void SMSSM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   m0 = pars(0);
   m12 = pars(1);
   TanBeta = pars(2);
   SignMu = pars(3);
   Azero = pars(4);
   LambdaInput = pars(5);
   KappaInput = pars(6);
   LambdaSInput = pars(7);
   L1Input = pars(8);
   MSInput = pars(9);
   BMSInput = pars(10);

}

std::ostream& operator<<(std::ostream& ostr, const SMSSM_input_parameters& input)
{
   ostr << "m0 = " << INPUT(m0) << ", ";
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "SignMu = " << INPUT(SignMu) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "LambdaInput = " << INPUT(LambdaInput) << ", ";
   ostr << "KappaInput = " << INPUT(KappaInput) << ", ";
   ostr << "LambdaSInput = " << INPUT(LambdaSInput) << ", ";
   ostr << "L1Input = " << INPUT(L1Input) << ", ";
   ostr << "MSInput = " << INPUT(MSInput) << ", ";
   ostr << "BMSInput = " << INPUT(BMSInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
