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

// File generated at Sun 26 Aug 2018 13:52:20

#include "CE6SSM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd CE6SSM_input_parameters::get() const
{
   Eigen::ArrayXd pars(10);

   pars(0) = TanBeta;
   pars(1) = m0SqGuess;
   pars(2) = m12Guess;
   pars(3) = AzeroGuess;
   pars(4) = LambdaInput;
   pars(5) = KappaInput;
   pars(6) = MuPrimeInput;
   pars(7) = BMuPrimeInput;
   pars(8) = vsInput;
   pars(9) = Lambda12Input;

   return pars;
}

void CE6SSM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   TanBeta = pars(0);
   m0SqGuess = pars(1);
   m12Guess = pars(2);
   AzeroGuess = pars(3);
   LambdaInput = pars(4);
   KappaInput = pars(5);
   MuPrimeInput = pars(6);
   BMuPrimeInput = pars(7);
   vsInput = pars(8);
   Lambda12Input = pars(9);

}

std::ostream& operator<<(std::ostream& ostr, const CE6SSM_input_parameters& input)
{
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "m0SqGuess = " << INPUT(m0SqGuess) << ", ";
   ostr << "m12Guess = " << INPUT(m12Guess) << ", ";
   ostr << "AzeroGuess = " << INPUT(AzeroGuess) << ", ";
   ostr << "LambdaInput = " << INPUT(LambdaInput) << ", ";
   ostr << "KappaInput = " << INPUT(KappaInput) << ", ";
   ostr << "MuPrimeInput = " << INPUT(MuPrimeInput) << ", ";
   ostr << "BMuPrimeInput = " << INPUT(BMuPrimeInput) << ", ";
   ostr << "vsInput = " << INPUT(vsInput) << ", ";
   ostr << "Lambda12Input = " << INPUT(Lambda12Input) << ", ";

   return ostr;
}

} // namespace flexiblesusy
