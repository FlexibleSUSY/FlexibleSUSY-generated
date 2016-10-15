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

// File generated at Sat 15 Oct 2016 15:34:28

#include "UMSSM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd UMSSM_input_parameters::get() const
{
   Eigen::ArrayXd pars(16);

   pars(0) = m0;
   pars(1) = m12;
   pars(2) = TanBeta;
   pars(3) = Azero;
   pars(4) = LambdaInput;
   pars(5) = ALambdaInput;
   pars(6) = vSInput;
   pars(7) = Qq;
   pars(8) = Ql;
   pars(9) = QHd;
   pars(10) = QHu;
   pars(11) = Qd;
   pars(12) = Qu;
   pars(13) = Qe;
   pars(14) = Qs;
   pars(15) = Qv;

   return pars;
}

void UMSSM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   m0 = pars(0);
   m12 = pars(1);
   TanBeta = pars(2);
   Azero = pars(3);
   LambdaInput = pars(4);
   ALambdaInput = pars(5);
   vSInput = pars(6);
   Qq = pars(7);
   Ql = pars(8);
   QHd = pars(9);
   QHu = pars(10);
   Qd = pars(11);
   Qu = pars(12);
   Qe = pars(13);
   Qs = pars(14);
   Qv = pars(15);

}

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
   ostr << "Qv = " << INPUT(Qv) << ", ";

   return ostr;
}

} // namespace flexiblesusy
