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

// File generated at Wed 29 Jun 2016 12:15:09

#include "E6SSM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd E6SSM_input_parameters::get() const
{
   Eigen::ArrayXd pars(10);

   pars(0) = m0;
   pars(1) = m12;
   pars(2) = TanBeta;
   pars(3) = Azero;
   pars(4) = LambdaInput;
   pars(5) = KappaInput;
   pars(6) = muPrimeInput;
   pars(7) = BmuPrimeInput;
   pars(8) = vSInput;
   pars(9) = Lambda12Input;

   return pars;
}

void E6SSM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   m0 = pars(0);
   m12 = pars(1);
   TanBeta = pars(2);
   Azero = pars(3);
   LambdaInput = pars(4);
   KappaInput = pars(5);
   muPrimeInput = pars(6);
   BmuPrimeInput = pars(7);
   vSInput = pars(8);
   Lambda12Input = pars(9);

}

std::ostream& operator<<(std::ostream& ostr, const E6SSM_input_parameters& input)
{
   ostr << "m0 = " << INPUT(m0) << ", ";
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "LambdaInput = " << INPUT(LambdaInput) << ", ";
   ostr << "KappaInput = " << INPUT(KappaInput) << ", ";
   ostr << "muPrimeInput = " << INPUT(muPrimeInput) << ", ";
   ostr << "BmuPrimeInput = " << INPUT(BmuPrimeInput) << ", ";
   ostr << "vSInput = " << INPUT(vSInput) << ", ";
   ostr << "Lambda12Input = " << INPUT(Lambda12Input) << ", ";

   return ostr;
}

} // namespace flexiblesusy
