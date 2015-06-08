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

// File generated at Mon 8 Jun 2015 17:56:28

#include "lowNMSSM_input_parameters.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

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
