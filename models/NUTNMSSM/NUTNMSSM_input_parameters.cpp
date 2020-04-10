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

// File generated at Fri 10 Apr 2020 20:27:06

#include "NUTNMSSM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd NUTNMSSM_input_parameters::get() const
{
   Eigen::ArrayXd pars(9);

   pars(0) = m0;
   pars(1) = m12;
   pars(2) = TanBeta;
   pars(3) = Azero;
   pars(4) = LambdaInput;
   pars(5) = KappaInput;
   pars(6) = ALambdaInput;
   pars(7) = AKappaInput;
   pars(8) = MuEff;

   return pars;
}

void NUTNMSSM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   m0 = pars(0);
   m12 = pars(1);
   TanBeta = pars(2);
   Azero = pars(3);
   LambdaInput = pars(4);
   KappaInput = pars(5);
   ALambdaInput = pars(6);
   AKappaInput = pars(7);
   MuEff = pars(8);

}

std::ostream& operator<<(std::ostream& ostr, const NUTNMSSM_input_parameters& input)
{
   ostr << "m0 = " << INPUT(m0) << ", ";
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "LambdaInput = " << INPUT(LambdaInput) << ", ";
   ostr << "KappaInput = " << INPUT(KappaInput) << ", ";
   ostr << "ALambdaInput = " << INPUT(ALambdaInput) << ", ";
   ostr << "AKappaInput = " << INPUT(AKappaInput) << ", ";
   ostr << "MuEff = " << INPUT(MuEff) << ", ";

   return ostr;
}

} // namespace flexiblesusy
