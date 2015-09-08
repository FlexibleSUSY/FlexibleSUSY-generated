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

// File generated at Tue 8 Sep 2015 12:13:34

#include "MRSSM_input_parameters.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

std::ostream& operator<<(std::ostream& ostr, const MRSSM_input_parameters& input)
{
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "LamTDInput = " << INPUT(LamTDInput) << ", ";
   ostr << "LamTUInput = " << INPUT(LamTUInput) << ", ";
   ostr << "LamSDInput = " << INPUT(LamSDInput) << ", ";
   ostr << "LamSUInput = " << INPUT(LamSUInput) << ", ";
   ostr << "MuInput = " << INPUT(MuInput) << ", ";
   ostr << "MuDInput = " << INPUT(MuDInput) << ", ";
   ostr << "MuUInput = " << INPUT(MuUInput) << ", ";
   ostr << "vTInput = " << INPUT(vTInput) << ", ";
   ostr << "vSInput = " << INPUT(vSInput) << ", ";
   ostr << "BMuInput = " << INPUT(BMuInput) << ", ";
   ostr << "BMuDInput = " << INPUT(BMuDInput) << ", ";
   ostr << "BMuUInput = " << INPUT(BMuUInput) << ", ";
   ostr << "mq2Input = " << INPUT(mq2Input) << ", ";
   ostr << "ml2Input = " << INPUT(ml2Input) << ", ";
   ostr << "md2Input = " << INPUT(md2Input) << ", ";
   ostr << "mu2Input = " << INPUT(mu2Input) << ", ";
   ostr << "me2Input = " << INPUT(me2Input) << ", ";
   ostr << "moc2Input = " << INPUT(moc2Input) << ", ";
   ostr << "mRd2Input = " << INPUT(mRd2Input) << ", ";
   ostr << "mRu2Input = " << INPUT(mRu2Input) << ", ";
   ostr << "MDBSInput = " << INPUT(MDBSInput) << ", ";
   ostr << "MDWBTInput = " << INPUT(MDWBTInput) << ", ";
   ostr << "MDGocInput = " << INPUT(MDGocInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
