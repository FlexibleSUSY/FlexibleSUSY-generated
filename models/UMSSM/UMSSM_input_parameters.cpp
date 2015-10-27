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

// File generated at Tue 27 Oct 2015 15:15:10

#include "UMSSM_input_parameters.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

std::ostream& operator<<(std::ostream& ostr, const UMSSM_input_parameters& input)
{
   ostr << "m0 = " << INPUT(m0) << ", ";
   ostr << "m12 = " << INPUT(m12) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "Azero = " << INPUT(Azero) << ", ";
   ostr << "LambdaInput = " << INPUT(LambdaInput) << ", ";
   ostr << "ALambdaInput = " << INPUT(ALambdaInput) << ", ";
   ostr << "vSInput = " << INPUT(vSInput) << ", ";
   ostr << "Qq = " << INPUT(Qq) << ", ";
   ostr << "Ql = " << INPUT(Ql) << ", ";
   ostr << "QHd = " << INPUT(QHd) << ", ";
   ostr << "QHu = " << INPUT(QHu) << ", ";
   ostr << "Qd = " << INPUT(Qd) << ", ";
   ostr << "Qu = " << INPUT(Qu) << ", ";
   ostr << "Qe = " << INPUT(Qe) << ", ";
   ostr << "Qs = " << INPUT(Qs) << ", ";

   return ostr;
}

} // namespace flexiblesusy
