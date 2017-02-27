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

// File generated at Mon 27 Feb 2017 13:35:48

#include "lowNMSSM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd lowNMSSM_input_parameters::get() const
{
   Eigen::ArrayXd pars(28);

   pars(0) = Qin;
   pars(1) = M1Input;
   pars(2) = M2Input;
   pars(3) = M3Input;
   pars(4) = AtInput;
   pars(5) = AbInput;
   pars(6) = ATauInput;
   pars(7) = TanBeta;
   pars(8) = ml1Input;
   pars(9) = ml2Input;
   pars(10) = ml3Input;
   pars(11) = me1Input;
   pars(12) = me2Input;
   pars(13) = me3Input;
   pars(14) = mq1Input;
   pars(15) = mq2Input;
   pars(16) = mq3Input;
   pars(17) = md1Input;
   pars(18) = md2Input;
   pars(19) = md3Input;
   pars(20) = mu1Input;
   pars(21) = mu2Input;
   pars(22) = mu3Input;
   pars(23) = LambdaInput;
   pars(24) = KappaInput;
   pars(25) = ALambdaInput;
   pars(26) = AKappaInput;
   pars(27) = MuEffInput;

   return pars;
}

void lowNMSSM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   Qin = pars(0);
   M1Input = pars(1);
   M2Input = pars(2);
   M3Input = pars(3);
   AtInput = pars(4);
   AbInput = pars(5);
   ATauInput = pars(6);
   TanBeta = pars(7);
   ml1Input = pars(8);
   ml2Input = pars(9);
   ml3Input = pars(10);
   me1Input = pars(11);
   me2Input = pars(12);
   me3Input = pars(13);
   mq1Input = pars(14);
   mq2Input = pars(15);
   mq3Input = pars(16);
   md1Input = pars(17);
   md2Input = pars(18);
   md3Input = pars(19);
   mu1Input = pars(20);
   mu2Input = pars(21);
   mu3Input = pars(22);
   LambdaInput = pars(23);
   KappaInput = pars(24);
   ALambdaInput = pars(25);
   AKappaInput = pars(26);
   MuEffInput = pars(27);

}

std::ostream& operator<<(std::ostream& ostr, const lowNMSSM_input_parameters& input)
{
   ostr << "Qin = " << INPUT(Qin) << ", ";
   ostr << "M1Input = " << INPUT(M1Input) << ", ";
   ostr << "M2Input = " << INPUT(M2Input) << ", ";
   ostr << "M3Input = " << INPUT(M3Input) << ", ";
   ostr << "AtInput = " << INPUT(AtInput) << ", ";
   ostr << "AbInput = " << INPUT(AbInput) << ", ";
   ostr << "ATauInput = " << INPUT(ATauInput) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "ml1Input = " << INPUT(ml1Input) << ", ";
   ostr << "ml2Input = " << INPUT(ml2Input) << ", ";
   ostr << "ml3Input = " << INPUT(ml3Input) << ", ";
   ostr << "me1Input = " << INPUT(me1Input) << ", ";
   ostr << "me2Input = " << INPUT(me2Input) << ", ";
   ostr << "me3Input = " << INPUT(me3Input) << ", ";
   ostr << "mq1Input = " << INPUT(mq1Input) << ", ";
   ostr << "mq2Input = " << INPUT(mq2Input) << ", ";
   ostr << "mq3Input = " << INPUT(mq3Input) << ", ";
   ostr << "md1Input = " << INPUT(md1Input) << ", ";
   ostr << "md2Input = " << INPUT(md2Input) << ", ";
   ostr << "md3Input = " << INPUT(md3Input) << ", ";
   ostr << "mu1Input = " << INPUT(mu1Input) << ", ";
   ostr << "mu2Input = " << INPUT(mu2Input) << ", ";
   ostr << "mu3Input = " << INPUT(mu3Input) << ", ";
   ostr << "LambdaInput = " << INPUT(LambdaInput) << ", ";
   ostr << "KappaInput = " << INPUT(KappaInput) << ", ";
   ostr << "ALambdaInput = " << INPUT(ALambdaInput) << ", ";
   ostr << "AKappaInput = " << INPUT(AKappaInput) << ", ";
   ostr << "MuEffInput = " << INPUT(MuEffInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
