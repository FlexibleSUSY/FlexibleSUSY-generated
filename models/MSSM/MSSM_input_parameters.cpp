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


#include "MSSM_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd MSSM_input_parameters::get() const
{
   Eigen::ArrayXd pars(80);

   pars(0) = TanBeta;
   pars(1) = SignMu;
   pars(2) = Qin;
   pars(3) = mHd2IN;
   pars(4) = mHu2IN;
   pars(5) = Aeij(0,0);
   pars(6) = Aeij(0,1);
   pars(7) = Aeij(0,2);
   pars(8) = Aeij(1,0);
   pars(9) = Aeij(1,1);
   pars(10) = Aeij(1,2);
   pars(11) = Aeij(2,0);
   pars(12) = Aeij(2,1);
   pars(13) = Aeij(2,2);
   pars(14) = Adij(0,0);
   pars(15) = Adij(0,1);
   pars(16) = Adij(0,2);
   pars(17) = Adij(1,0);
   pars(18) = Adij(1,1);
   pars(19) = Adij(1,2);
   pars(20) = Adij(2,0);
   pars(21) = Adij(2,1);
   pars(22) = Adij(2,2);
   pars(23) = Auij(0,0);
   pars(24) = Auij(0,1);
   pars(25) = Auij(0,2);
   pars(26) = Auij(1,0);
   pars(27) = Auij(1,1);
   pars(28) = Auij(1,2);
   pars(29) = Auij(2,0);
   pars(30) = Auij(2,1);
   pars(31) = Auij(2,2);
   pars(32) = mq2Input(0,0);
   pars(33) = mq2Input(0,1);
   pars(34) = mq2Input(0,2);
   pars(35) = mq2Input(1,0);
   pars(36) = mq2Input(1,1);
   pars(37) = mq2Input(1,2);
   pars(38) = mq2Input(2,0);
   pars(39) = mq2Input(2,1);
   pars(40) = mq2Input(2,2);
   pars(41) = ml2Input(0,0);
   pars(42) = ml2Input(0,1);
   pars(43) = ml2Input(0,2);
   pars(44) = ml2Input(1,0);
   pars(45) = ml2Input(1,1);
   pars(46) = ml2Input(1,2);
   pars(47) = ml2Input(2,0);
   pars(48) = ml2Input(2,1);
   pars(49) = ml2Input(2,2);
   pars(50) = md2Input(0,0);
   pars(51) = md2Input(0,1);
   pars(52) = md2Input(0,2);
   pars(53) = md2Input(1,0);
   pars(54) = md2Input(1,1);
   pars(55) = md2Input(1,2);
   pars(56) = md2Input(2,0);
   pars(57) = md2Input(2,1);
   pars(58) = md2Input(2,2);
   pars(59) = mu2Input(0,0);
   pars(60) = mu2Input(0,1);
   pars(61) = mu2Input(0,2);
   pars(62) = mu2Input(1,0);
   pars(63) = mu2Input(1,1);
   pars(64) = mu2Input(1,2);
   pars(65) = mu2Input(2,0);
   pars(66) = mu2Input(2,1);
   pars(67) = mu2Input(2,2);
   pars(68) = me2Input(0,0);
   pars(69) = me2Input(0,1);
   pars(70) = me2Input(0,2);
   pars(71) = me2Input(1,0);
   pars(72) = me2Input(1,1);
   pars(73) = me2Input(1,2);
   pars(74) = me2Input(2,0);
   pars(75) = me2Input(2,1);
   pars(76) = me2Input(2,2);
   pars(77) = MassBInput;
   pars(78) = MassWBInput;
   pars(79) = MassGInput;

   return pars;
}

void MSSM_input_parameters::set(const Eigen::ArrayXd& pars)
{
   TanBeta = pars(0);
   SignMu = pars(1);
   Qin = pars(2);
   mHd2IN = pars(3);
   mHu2IN = pars(4);
   Aeij(0,0) = pars(5);
   Aeij(0,1) = pars(6);
   Aeij(0,2) = pars(7);
   Aeij(1,0) = pars(8);
   Aeij(1,1) = pars(9);
   Aeij(1,2) = pars(10);
   Aeij(2,0) = pars(11);
   Aeij(2,1) = pars(12);
   Aeij(2,2) = pars(13);
   Adij(0,0) = pars(14);
   Adij(0,1) = pars(15);
   Adij(0,2) = pars(16);
   Adij(1,0) = pars(17);
   Adij(1,1) = pars(18);
   Adij(1,2) = pars(19);
   Adij(2,0) = pars(20);
   Adij(2,1) = pars(21);
   Adij(2,2) = pars(22);
   Auij(0,0) = pars(23);
   Auij(0,1) = pars(24);
   Auij(0,2) = pars(25);
   Auij(1,0) = pars(26);
   Auij(1,1) = pars(27);
   Auij(1,2) = pars(28);
   Auij(2,0) = pars(29);
   Auij(2,1) = pars(30);
   Auij(2,2) = pars(31);
   mq2Input(0,0) = pars(32);
   mq2Input(0,1) = pars(33);
   mq2Input(0,2) = pars(34);
   mq2Input(1,0) = pars(35);
   mq2Input(1,1) = pars(36);
   mq2Input(1,2) = pars(37);
   mq2Input(2,0) = pars(38);
   mq2Input(2,1) = pars(39);
   mq2Input(2,2) = pars(40);
   ml2Input(0,0) = pars(41);
   ml2Input(0,1) = pars(42);
   ml2Input(0,2) = pars(43);
   ml2Input(1,0) = pars(44);
   ml2Input(1,1) = pars(45);
   ml2Input(1,2) = pars(46);
   ml2Input(2,0) = pars(47);
   ml2Input(2,1) = pars(48);
   ml2Input(2,2) = pars(49);
   md2Input(0,0) = pars(50);
   md2Input(0,1) = pars(51);
   md2Input(0,2) = pars(52);
   md2Input(1,0) = pars(53);
   md2Input(1,1) = pars(54);
   md2Input(1,2) = pars(55);
   md2Input(2,0) = pars(56);
   md2Input(2,1) = pars(57);
   md2Input(2,2) = pars(58);
   mu2Input(0,0) = pars(59);
   mu2Input(0,1) = pars(60);
   mu2Input(0,2) = pars(61);
   mu2Input(1,0) = pars(62);
   mu2Input(1,1) = pars(63);
   mu2Input(1,2) = pars(64);
   mu2Input(2,0) = pars(65);
   mu2Input(2,1) = pars(66);
   mu2Input(2,2) = pars(67);
   me2Input(0,0) = pars(68);
   me2Input(0,1) = pars(69);
   me2Input(0,2) = pars(70);
   me2Input(1,0) = pars(71);
   me2Input(1,1) = pars(72);
   me2Input(1,2) = pars(73);
   me2Input(2,0) = pars(74);
   me2Input(2,1) = pars(75);
   me2Input(2,2) = pars(76);
   MassBInput = pars(77);
   MassWBInput = pars(78);
   MassGInput = pars(79);

}

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
