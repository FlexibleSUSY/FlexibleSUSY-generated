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

// File generated at Thu 15 Dec 2016 12:39:44

#include "MRSSMtower_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd MRSSMtower_input_parameters::get() const
{
   Eigen::ArrayXd pars(62);

   pars(0) = TanBeta;
   pars(1) = MS;
   pars(2) = LamTDInput;
   pars(3) = LamTUInput;
   pars(4) = LamSDInput;
   pars(5) = LamSUInput;
   pars(6) = MuDInput;
   pars(7) = MuUInput;
   pars(8) = BMuInput;
   pars(9) = mq2Input(0,0);
   pars(10) = mq2Input(0,1);
   pars(11) = mq2Input(0,2);
   pars(12) = mq2Input(1,0);
   pars(13) = mq2Input(1,1);
   pars(14) = mq2Input(1,2);
   pars(15) = mq2Input(2,0);
   pars(16) = mq2Input(2,1);
   pars(17) = mq2Input(2,2);
   pars(18) = ml2Input(0,0);
   pars(19) = ml2Input(0,1);
   pars(20) = ml2Input(0,2);
   pars(21) = ml2Input(1,0);
   pars(22) = ml2Input(1,1);
   pars(23) = ml2Input(1,2);
   pars(24) = ml2Input(2,0);
   pars(25) = ml2Input(2,1);
   pars(26) = ml2Input(2,2);
   pars(27) = md2Input(0,0);
   pars(28) = md2Input(0,1);
   pars(29) = md2Input(0,2);
   pars(30) = md2Input(1,0);
   pars(31) = md2Input(1,1);
   pars(32) = md2Input(1,2);
   pars(33) = md2Input(2,0);
   pars(34) = md2Input(2,1);
   pars(35) = md2Input(2,2);
   pars(36) = mu2Input(0,0);
   pars(37) = mu2Input(0,1);
   pars(38) = mu2Input(0,2);
   pars(39) = mu2Input(1,0);
   pars(40) = mu2Input(1,1);
   pars(41) = mu2Input(1,2);
   pars(42) = mu2Input(2,0);
   pars(43) = mu2Input(2,1);
   pars(44) = mu2Input(2,2);
   pars(45) = me2Input(0,0);
   pars(46) = me2Input(0,1);
   pars(47) = me2Input(0,2);
   pars(48) = me2Input(1,0);
   pars(49) = me2Input(1,1);
   pars(50) = me2Input(1,2);
   pars(51) = me2Input(2,0);
   pars(52) = me2Input(2,1);
   pars(53) = me2Input(2,2);
   pars(54) = mS2Input;
   pars(55) = mT2Input;
   pars(56) = moc2Input;
   pars(57) = mRd2Input;
   pars(58) = mRu2Input;
   pars(59) = MDBSInput;
   pars(60) = MDWBTInput;
   pars(61) = MDGocInput;

   return pars;
}

void MRSSMtower_input_parameters::set(const Eigen::ArrayXd& pars)
{
   TanBeta = pars(0);
   MS = pars(1);
   LamTDInput = pars(2);
   LamTUInput = pars(3);
   LamSDInput = pars(4);
   LamSUInput = pars(5);
   MuDInput = pars(6);
   MuUInput = pars(7);
   BMuInput = pars(8);
   mq2Input(0,0) = pars(9);
   mq2Input(0,1) = pars(10);
   mq2Input(0,2) = pars(11);
   mq2Input(1,0) = pars(12);
   mq2Input(1,1) = pars(13);
   mq2Input(1,2) = pars(14);
   mq2Input(2,0) = pars(15);
   mq2Input(2,1) = pars(16);
   mq2Input(2,2) = pars(17);
   ml2Input(0,0) = pars(18);
   ml2Input(0,1) = pars(19);
   ml2Input(0,2) = pars(20);
   ml2Input(1,0) = pars(21);
   ml2Input(1,1) = pars(22);
   ml2Input(1,2) = pars(23);
   ml2Input(2,0) = pars(24);
   ml2Input(2,1) = pars(25);
   ml2Input(2,2) = pars(26);
   md2Input(0,0) = pars(27);
   md2Input(0,1) = pars(28);
   md2Input(0,2) = pars(29);
   md2Input(1,0) = pars(30);
   md2Input(1,1) = pars(31);
   md2Input(1,2) = pars(32);
   md2Input(2,0) = pars(33);
   md2Input(2,1) = pars(34);
   md2Input(2,2) = pars(35);
   mu2Input(0,0) = pars(36);
   mu2Input(0,1) = pars(37);
   mu2Input(0,2) = pars(38);
   mu2Input(1,0) = pars(39);
   mu2Input(1,1) = pars(40);
   mu2Input(1,2) = pars(41);
   mu2Input(2,0) = pars(42);
   mu2Input(2,1) = pars(43);
   mu2Input(2,2) = pars(44);
   me2Input(0,0) = pars(45);
   me2Input(0,1) = pars(46);
   me2Input(0,2) = pars(47);
   me2Input(1,0) = pars(48);
   me2Input(1,1) = pars(49);
   me2Input(1,2) = pars(50);
   me2Input(2,0) = pars(51);
   me2Input(2,1) = pars(52);
   me2Input(2,2) = pars(53);
   mS2Input = pars(54);
   mT2Input = pars(55);
   moc2Input = pars(56);
   mRd2Input = pars(57);
   mRu2Input = pars(58);
   MDBSInput = pars(59);
   MDWBTInput = pars(60);
   MDGocInput = pars(61);

}

std::ostream& operator<<(std::ostream& ostr, const MRSSMtower_input_parameters& input)
{
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "MS = " << INPUT(MS) << ", ";
   ostr << "LamTDInput = " << INPUT(LamTDInput) << ", ";
   ostr << "LamTUInput = " << INPUT(LamTUInput) << ", ";
   ostr << "LamSDInput = " << INPUT(LamSDInput) << ", ";
   ostr << "LamSUInput = " << INPUT(LamSUInput) << ", ";
   ostr << "MuDInput = " << INPUT(MuDInput) << ", ";
   ostr << "MuUInput = " << INPUT(MuUInput) << ", ";
   ostr << "BMuInput = " << INPUT(BMuInput) << ", ";
   ostr << "mq2Input = " << INPUT(mq2Input) << ", ";
   ostr << "ml2Input = " << INPUT(ml2Input) << ", ";
   ostr << "md2Input = " << INPUT(md2Input) << ", ";
   ostr << "mu2Input = " << INPUT(mu2Input) << ", ";
   ostr << "me2Input = " << INPUT(me2Input) << ", ";
   ostr << "mS2Input = " << INPUT(mS2Input) << ", ";
   ostr << "mT2Input = " << INPUT(mT2Input) << ", ";
   ostr << "moc2Input = " << INPUT(moc2Input) << ", ";
   ostr << "mRd2Input = " << INPUT(mRd2Input) << ", ";
   ostr << "mRu2Input = " << INPUT(mRu2Input) << ", ";
   ostr << "MDBSInput = " << INPUT(MDBSInput) << ", ";
   ostr << "MDWBTInput = " << INPUT(MDWBTInput) << ", ";
   ostr << "MDGocInput = " << INPUT(MDGocInput) << ", ";

   return ostr;
}

} // namespace flexiblesusy
