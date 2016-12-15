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

// File generated at Thu 15 Dec 2016 12:46:25

#include "MRSSM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd MRSSM_input_parameters::get() const
{
   Eigen::ArrayXd pars(64);

   pars(0) = TanBeta;
   pars(1) = LamTDInput;
   pars(2) = LamTUInput;
   pars(3) = LamSDInput;
   pars(4) = LamSUInput;
   pars(5) = MuInput;
   pars(6) = MuDInput;
   pars(7) = MuUInput;
   pars(8) = vTInput;
   pars(9) = vSInput;
   pars(10) = BMuInput;
   pars(11) = BMuDInput;
   pars(12) = BMuUInput;
   pars(13) = mq2Input(0,0);
   pars(14) = mq2Input(0,1);
   pars(15) = mq2Input(0,2);
   pars(16) = mq2Input(1,0);
   pars(17) = mq2Input(1,1);
   pars(18) = mq2Input(1,2);
   pars(19) = mq2Input(2,0);
   pars(20) = mq2Input(2,1);
   pars(21) = mq2Input(2,2);
   pars(22) = ml2Input(0,0);
   pars(23) = ml2Input(0,1);
   pars(24) = ml2Input(0,2);
   pars(25) = ml2Input(1,0);
   pars(26) = ml2Input(1,1);
   pars(27) = ml2Input(1,2);
   pars(28) = ml2Input(2,0);
   pars(29) = ml2Input(2,1);
   pars(30) = ml2Input(2,2);
   pars(31) = md2Input(0,0);
   pars(32) = md2Input(0,1);
   pars(33) = md2Input(0,2);
   pars(34) = md2Input(1,0);
   pars(35) = md2Input(1,1);
   pars(36) = md2Input(1,2);
   pars(37) = md2Input(2,0);
   pars(38) = md2Input(2,1);
   pars(39) = md2Input(2,2);
   pars(40) = mu2Input(0,0);
   pars(41) = mu2Input(0,1);
   pars(42) = mu2Input(0,2);
   pars(43) = mu2Input(1,0);
   pars(44) = mu2Input(1,1);
   pars(45) = mu2Input(1,2);
   pars(46) = mu2Input(2,0);
   pars(47) = mu2Input(2,1);
   pars(48) = mu2Input(2,2);
   pars(49) = me2Input(0,0);
   pars(50) = me2Input(0,1);
   pars(51) = me2Input(0,2);
   pars(52) = me2Input(1,0);
   pars(53) = me2Input(1,1);
   pars(54) = me2Input(1,2);
   pars(55) = me2Input(2,0);
   pars(56) = me2Input(2,1);
   pars(57) = me2Input(2,2);
   pars(58) = moc2Input;
   pars(59) = mRd2Input;
   pars(60) = mRu2Input;
   pars(61) = MDBSInput;
   pars(62) = MDWBTInput;
   pars(63) = MDGocInput;

   return pars;
}

void MRSSM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   TanBeta = pars(0);
   LamTDInput = pars(1);
   LamTUInput = pars(2);
   LamSDInput = pars(3);
   LamSUInput = pars(4);
   MuInput = pars(5);
   MuDInput = pars(6);
   MuUInput = pars(7);
   vTInput = pars(8);
   vSInput = pars(9);
   BMuInput = pars(10);
   BMuDInput = pars(11);
   BMuUInput = pars(12);
   mq2Input(0,0) = pars(13);
   mq2Input(0,1) = pars(14);
   mq2Input(0,2) = pars(15);
   mq2Input(1,0) = pars(16);
   mq2Input(1,1) = pars(17);
   mq2Input(1,2) = pars(18);
   mq2Input(2,0) = pars(19);
   mq2Input(2,1) = pars(20);
   mq2Input(2,2) = pars(21);
   ml2Input(0,0) = pars(22);
   ml2Input(0,1) = pars(23);
   ml2Input(0,2) = pars(24);
   ml2Input(1,0) = pars(25);
   ml2Input(1,1) = pars(26);
   ml2Input(1,2) = pars(27);
   ml2Input(2,0) = pars(28);
   ml2Input(2,1) = pars(29);
   ml2Input(2,2) = pars(30);
   md2Input(0,0) = pars(31);
   md2Input(0,1) = pars(32);
   md2Input(0,2) = pars(33);
   md2Input(1,0) = pars(34);
   md2Input(1,1) = pars(35);
   md2Input(1,2) = pars(36);
   md2Input(2,0) = pars(37);
   md2Input(2,1) = pars(38);
   md2Input(2,2) = pars(39);
   mu2Input(0,0) = pars(40);
   mu2Input(0,1) = pars(41);
   mu2Input(0,2) = pars(42);
   mu2Input(1,0) = pars(43);
   mu2Input(1,1) = pars(44);
   mu2Input(1,2) = pars(45);
   mu2Input(2,0) = pars(46);
   mu2Input(2,1) = pars(47);
   mu2Input(2,2) = pars(48);
   me2Input(0,0) = pars(49);
   me2Input(0,1) = pars(50);
   me2Input(0,2) = pars(51);
   me2Input(1,0) = pars(52);
   me2Input(1,1) = pars(53);
   me2Input(1,2) = pars(54);
   me2Input(2,0) = pars(55);
   me2Input(2,1) = pars(56);
   me2Input(2,2) = pars(57);
   moc2Input = pars(58);
   mRd2Input = pars(59);
   mRu2Input = pars(60);
   MDBSInput = pars(61);
   MDWBTInput = pars(62);
   MDGocInput = pars(63);

}

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
