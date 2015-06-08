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

// File generated at Mon 8 Jun 2015 18:08:43

#include "MSSM_input_parameters.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

std::ostream& operator<<(std::ostream& ostr, const MSSM_input_parameters& input)
{
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "SignMu = " << INPUT(SignMu) << ", ";
   ostr << "Qin = " << INPUT(Qin) << ", ";
   ostr << "mHd2IN = " << INPUT(mHd2IN) << ", ";
   ostr << "mHu2IN = " << INPUT(mHu2IN) << ", ";
   ostr << "Aeij = " << INPUT(Aeij) << ", ";
   ostr << "Adij = " << INPUT(Adij) << ", ";
   ostr << "Auij = " << INPUT(Auij) << ", ";
   ostr << "mq2Input = " << INPUT(mq2Input) << ", ";
   ostr << "ml2Input = " << INPUT(ml2Input) << ", ";
   ostr << "md2Input = " << INPUT(md2Input) << ", ";
   ostr << "mu2Input = " << INPUT(mu2Input) << ", ";
   ostr << "me2Input = " << INPUT(me2Input) << ", ";
   ostr << "MassBInput = " << INPUT(MassBInput) << ", ";
   ostr << "MassWBInput = " << INPUT(MassWBInput) << ", ";
   ostr << "MassGInput = " << INPUT(MassGInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
