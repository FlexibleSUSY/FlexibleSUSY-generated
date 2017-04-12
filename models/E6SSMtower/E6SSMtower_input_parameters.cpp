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

// File generated at Wed 12 Apr 2017 11:27:18

#include "E6SSMtower_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd E6SSMtower_input_parameters::get() const
{
   Eigen::ArrayXd pars(142);

   pars(0) = MSUSY;
   pars(1) = M1Input;
   pars(2) = M2Input;
   pars(3) = M3Input;
   pars(4) = MuInput;
   pars(5) = mAInput;
   pars(6) = TanBeta;
   pars(7) = LambdaInput;
   pars(8) = gNInput;
   pars(9) = M4Input;
   pars(10) = mHp2Input;
   pars(11) = mHpbar2Input;
   pars(12) = MuPrInput;
   pars(13) = BMuPrInput;
   pars(14) = mq2Input(0,0);
   pars(15) = mq2Input(0,1);
   pars(16) = mq2Input(0,2);
   pars(17) = mq2Input(1,0);
   pars(18) = mq2Input(1,1);
   pars(19) = mq2Input(1,2);
   pars(20) = mq2Input(2,0);
   pars(21) = mq2Input(2,1);
   pars(22) = mq2Input(2,2);
   pars(23) = mu2Input(0,0);
   pars(24) = mu2Input(0,1);
   pars(25) = mu2Input(0,2);
   pars(26) = mu2Input(1,0);
   pars(27) = mu2Input(1,1);
   pars(28) = mu2Input(1,2);
   pars(29) = mu2Input(2,0);
   pars(30) = mu2Input(2,1);
   pars(31) = mu2Input(2,2);
   pars(32) = md2Input(0,0);
   pars(33) = md2Input(0,1);
   pars(34) = md2Input(0,2);
   pars(35) = md2Input(1,0);
   pars(36) = md2Input(1,1);
   pars(37) = md2Input(1,2);
   pars(38) = md2Input(2,0);
   pars(39) = md2Input(2,1);
   pars(40) = md2Input(2,2);
   pars(41) = ml2Input(0,0);
   pars(42) = ml2Input(0,1);
   pars(43) = ml2Input(0,2);
   pars(44) = ml2Input(1,0);
   pars(45) = ml2Input(1,1);
   pars(46) = ml2Input(1,2);
   pars(47) = ml2Input(2,0);
   pars(48) = ml2Input(2,1);
   pars(49) = ml2Input(2,2);
   pars(50) = me2Input(0,0);
   pars(51) = me2Input(0,1);
   pars(52) = me2Input(0,2);
   pars(53) = me2Input(1,0);
   pars(54) = me2Input(1,1);
   pars(55) = me2Input(1,2);
   pars(56) = me2Input(2,0);
   pars(57) = me2Input(2,1);
   pars(58) = me2Input(2,2);
   pars(59) = AuInput(0,0);
   pars(60) = AuInput(0,1);
   pars(61) = AuInput(0,2);
   pars(62) = AuInput(1,0);
   pars(63) = AuInput(1,1);
   pars(64) = AuInput(1,2);
   pars(65) = AuInput(2,0);
   pars(66) = AuInput(2,1);
   pars(67) = AuInput(2,2);
   pars(68) = AdInput(0,0);
   pars(69) = AdInput(0,1);
   pars(70) = AdInput(0,2);
   pars(71) = AdInput(1,0);
   pars(72) = AdInput(1,1);
   pars(73) = AdInput(1,2);
   pars(74) = AdInput(2,0);
   pars(75) = AdInput(2,1);
   pars(76) = AdInput(2,2);
   pars(77) = AeInput(0,0);
   pars(78) = AeInput(0,1);
   pars(79) = AeInput(0,2);
   pars(80) = AeInput(1,0);
   pars(81) = AeInput(1,1);
   pars(82) = AeInput(1,2);
   pars(83) = AeInput(2,0);
   pars(84) = AeInput(2,1);
   pars(85) = AeInput(2,2);
   pars(86) = Lambda12Input(0,0);
   pars(87) = Lambda12Input(0,1);
   pars(88) = Lambda12Input(1,0);
   pars(89) = Lambda12Input(1,1);
   pars(90) = ALambda12Input(0,0);
   pars(91) = ALambda12Input(0,1);
   pars(92) = ALambda12Input(1,0);
   pars(93) = ALambda12Input(1,1);
   pars(94) = KappaInput(0,0);
   pars(95) = KappaInput(0,1);
   pars(96) = KappaInput(0,2);
   pars(97) = KappaInput(1,0);
   pars(98) = KappaInput(1,1);
   pars(99) = KappaInput(1,2);
   pars(100) = KappaInput(2,0);
   pars(101) = KappaInput(2,1);
   pars(102) = KappaInput(2,2);
   pars(103) = AKappaInput(0,0);
   pars(104) = AKappaInput(0,1);
   pars(105) = AKappaInput(0,2);
   pars(106) = AKappaInput(1,0);
   pars(107) = AKappaInput(1,1);
   pars(108) = AKappaInput(1,2);
   pars(109) = AKappaInput(2,0);
   pars(110) = AKappaInput(2,1);
   pars(111) = AKappaInput(2,2);
   pars(112) = mDx2Input(0,0);
   pars(113) = mDx2Input(0,1);
   pars(114) = mDx2Input(0,2);
   pars(115) = mDx2Input(1,0);
   pars(116) = mDx2Input(1,1);
   pars(117) = mDx2Input(1,2);
   pars(118) = mDx2Input(2,0);
   pars(119) = mDx2Input(2,1);
   pars(120) = mDx2Input(2,2);
   pars(121) = mDxbar2Input(0,0);
   pars(122) = mDxbar2Input(0,1);
   pars(123) = mDxbar2Input(0,2);
   pars(124) = mDxbar2Input(1,0);
   pars(125) = mDxbar2Input(1,1);
   pars(126) = mDxbar2Input(1,2);
   pars(127) = mDxbar2Input(2,0);
   pars(128) = mDxbar2Input(2,1);
   pars(129) = mDxbar2Input(2,2);
   pars(130) = mH1I2Input(0,0);
   pars(131) = mH1I2Input(0,1);
   pars(132) = mH1I2Input(1,0);
   pars(133) = mH1I2Input(1,1);
   pars(134) = mH2I2Input(0,0);
   pars(135) = mH2I2Input(0,1);
   pars(136) = mH2I2Input(1,0);
   pars(137) = mH2I2Input(1,1);
   pars(138) = msI2Input(0,0);
   pars(139) = msI2Input(0,1);
   pars(140) = msI2Input(1,0);
   pars(141) = msI2Input(1,1);

   return pars;
}

void E6SSMtower_input_parameters::set(const Eigen::ArrayXd& pars)
{
   MSUSY = pars(0);
   M1Input = pars(1);
   M2Input = pars(2);
   M3Input = pars(3);
   MuInput = pars(4);
   mAInput = pars(5);
   TanBeta = pars(6);
   LambdaInput = pars(7);
   gNInput = pars(8);
   M4Input = pars(9);
   mHp2Input = pars(10);
   mHpbar2Input = pars(11);
   MuPrInput = pars(12);
   BMuPrInput = pars(13);
   mq2Input(0,0) = pars(14);
   mq2Input(0,1) = pars(15);
   mq2Input(0,2) = pars(16);
   mq2Input(1,0) = pars(17);
   mq2Input(1,1) = pars(18);
   mq2Input(1,2) = pars(19);
   mq2Input(2,0) = pars(20);
   mq2Input(2,1) = pars(21);
   mq2Input(2,2) = pars(22);
   mu2Input(0,0) = pars(23);
   mu2Input(0,1) = pars(24);
   mu2Input(0,2) = pars(25);
   mu2Input(1,0) = pars(26);
   mu2Input(1,1) = pars(27);
   mu2Input(1,2) = pars(28);
   mu2Input(2,0) = pars(29);
   mu2Input(2,1) = pars(30);
   mu2Input(2,2) = pars(31);
   md2Input(0,0) = pars(32);
   md2Input(0,1) = pars(33);
   md2Input(0,2) = pars(34);
   md2Input(1,0) = pars(35);
   md2Input(1,1) = pars(36);
   md2Input(1,2) = pars(37);
   md2Input(2,0) = pars(38);
   md2Input(2,1) = pars(39);
   md2Input(2,2) = pars(40);
   ml2Input(0,0) = pars(41);
   ml2Input(0,1) = pars(42);
   ml2Input(0,2) = pars(43);
   ml2Input(1,0) = pars(44);
   ml2Input(1,1) = pars(45);
   ml2Input(1,2) = pars(46);
   ml2Input(2,0) = pars(47);
   ml2Input(2,1) = pars(48);
   ml2Input(2,2) = pars(49);
   me2Input(0,0) = pars(50);
   me2Input(0,1) = pars(51);
   me2Input(0,2) = pars(52);
   me2Input(1,0) = pars(53);
   me2Input(1,1) = pars(54);
   me2Input(1,2) = pars(55);
   me2Input(2,0) = pars(56);
   me2Input(2,1) = pars(57);
   me2Input(2,2) = pars(58);
   AuInput(0,0) = pars(59);
   AuInput(0,1) = pars(60);
   AuInput(0,2) = pars(61);
   AuInput(1,0) = pars(62);
   AuInput(1,1) = pars(63);
   AuInput(1,2) = pars(64);
   AuInput(2,0) = pars(65);
   AuInput(2,1) = pars(66);
   AuInput(2,2) = pars(67);
   AdInput(0,0) = pars(68);
   AdInput(0,1) = pars(69);
   AdInput(0,2) = pars(70);
   AdInput(1,0) = pars(71);
   AdInput(1,1) = pars(72);
   AdInput(1,2) = pars(73);
   AdInput(2,0) = pars(74);
   AdInput(2,1) = pars(75);
   AdInput(2,2) = pars(76);
   AeInput(0,0) = pars(77);
   AeInput(0,1) = pars(78);
   AeInput(0,2) = pars(79);
   AeInput(1,0) = pars(80);
   AeInput(1,1) = pars(81);
   AeInput(1,2) = pars(82);
   AeInput(2,0) = pars(83);
   AeInput(2,1) = pars(84);
   AeInput(2,2) = pars(85);
   Lambda12Input(0,0) = pars(86);
   Lambda12Input(0,1) = pars(87);
   Lambda12Input(1,0) = pars(88);
   Lambda12Input(1,1) = pars(89);
   ALambda12Input(0,0) = pars(90);
   ALambda12Input(0,1) = pars(91);
   ALambda12Input(1,0) = pars(92);
   ALambda12Input(1,1) = pars(93);
   KappaInput(0,0) = pars(94);
   KappaInput(0,1) = pars(95);
   KappaInput(0,2) = pars(96);
   KappaInput(1,0) = pars(97);
   KappaInput(1,1) = pars(98);
   KappaInput(1,2) = pars(99);
   KappaInput(2,0) = pars(100);
   KappaInput(2,1) = pars(101);
   KappaInput(2,2) = pars(102);
   AKappaInput(0,0) = pars(103);
   AKappaInput(0,1) = pars(104);
   AKappaInput(0,2) = pars(105);
   AKappaInput(1,0) = pars(106);
   AKappaInput(1,1) = pars(107);
   AKappaInput(1,2) = pars(108);
   AKappaInput(2,0) = pars(109);
   AKappaInput(2,1) = pars(110);
   AKappaInput(2,2) = pars(111);
   mDx2Input(0,0) = pars(112);
   mDx2Input(0,1) = pars(113);
   mDx2Input(0,2) = pars(114);
   mDx2Input(1,0) = pars(115);
   mDx2Input(1,1) = pars(116);
   mDx2Input(1,2) = pars(117);
   mDx2Input(2,0) = pars(118);
   mDx2Input(2,1) = pars(119);
   mDx2Input(2,2) = pars(120);
   mDxbar2Input(0,0) = pars(121);
   mDxbar2Input(0,1) = pars(122);
   mDxbar2Input(0,2) = pars(123);
   mDxbar2Input(1,0) = pars(124);
   mDxbar2Input(1,1) = pars(125);
   mDxbar2Input(1,2) = pars(126);
   mDxbar2Input(2,0) = pars(127);
   mDxbar2Input(2,1) = pars(128);
   mDxbar2Input(2,2) = pars(129);
   mH1I2Input(0,0) = pars(130);
   mH1I2Input(0,1) = pars(131);
   mH1I2Input(1,0) = pars(132);
   mH1I2Input(1,1) = pars(133);
   mH2I2Input(0,0) = pars(134);
   mH2I2Input(0,1) = pars(135);
   mH2I2Input(1,0) = pars(136);
   mH2I2Input(1,1) = pars(137);
   msI2Input(0,0) = pars(138);
   msI2Input(0,1) = pars(139);
   msI2Input(1,0) = pars(140);
   msI2Input(1,1) = pars(141);

}

std::ostream& operator<<(std::ostream& ostr, const E6SSMtower_input_parameters& input)
{
   ostr << "MSUSY = " << INPUT(MSUSY) << ", ";
   ostr << "M1Input = " << INPUT(M1Input) << ", ";
   ostr << "M2Input = " << INPUT(M2Input) << ", ";
   ostr << "M3Input = " << INPUT(M3Input) << ", ";
   ostr << "MuInput = " << INPUT(MuInput) << ", ";
   ostr << "mAInput = " << INPUT(mAInput) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "LambdaInput = " << INPUT(LambdaInput) << ", ";
   ostr << "gNInput = " << INPUT(gNInput) << ", ";
   ostr << "M4Input = " << INPUT(M4Input) << ", ";
   ostr << "mHp2Input = " << INPUT(mHp2Input) << ", ";
   ostr << "mHpbar2Input = " << INPUT(mHpbar2Input) << ", ";
   ostr << "MuPrInput = " << INPUT(MuPrInput) << ", ";
   ostr << "BMuPrInput = " << INPUT(BMuPrInput) << ", ";
   ostr << "mq2Input = " << INPUT(mq2Input) << ", ";
   ostr << "mu2Input = " << INPUT(mu2Input) << ", ";
   ostr << "md2Input = " << INPUT(md2Input) << ", ";
   ostr << "ml2Input = " << INPUT(ml2Input) << ", ";
   ostr << "me2Input = " << INPUT(me2Input) << ", ";
   ostr << "AuInput = " << INPUT(AuInput) << ", ";
   ostr << "AdInput = " << INPUT(AdInput) << ", ";
   ostr << "AeInput = " << INPUT(AeInput) << ", ";
   ostr << "Lambda12Input = " << INPUT(Lambda12Input) << ", ";
   ostr << "ALambda12Input = " << INPUT(ALambda12Input) << ", ";
   ostr << "KappaInput = " << INPUT(KappaInput) << ", ";
   ostr << "AKappaInput = " << INPUT(AKappaInput) << ", ";
   ostr << "mDx2Input = " << INPUT(mDx2Input) << ", ";
   ostr << "mDxbar2Input = " << INPUT(mDxbar2Input) << ", ";
   ostr << "mH1I2Input = " << INPUT(mH1I2Input) << ", ";
   ostr << "mH2I2Input = " << INPUT(mH2I2Input) << ", ";
   ostr << "msI2Input = " << INPUT(msI2Input) << ", ";

   return ostr;
}

} // namespace flexiblesusy
