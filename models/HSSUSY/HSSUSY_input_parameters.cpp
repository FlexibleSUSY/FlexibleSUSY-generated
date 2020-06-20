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


#include "HSSUSY_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd HSSUSY_input_parameters::get() const
{
   Eigen::ArrayXd pars(68);

   pars(0) = MSUSY;
   pars(1) = M1Input;
   pars(2) = M2Input;
   pars(3) = M3Input;
   pars(4) = MuInput;
   pars(5) = mAInput;
   pars(6) = MEWSB;
   pars(7) = AtInput;
   pars(8) = AbInput;
   pars(9) = AtauInput;
   pars(10) = TanBeta;
   pars(11) = LambdaLoopOrder;
   pars(12) = TwoLoopAtAs;
   pars(13) = TwoLoopAbAs;
   pars(14) = TwoLoopAtAb;
   pars(15) = TwoLoopAtauAtau;
   pars(16) = TwoLoopAtAt;
   pars(17) = DeltaEFT;
   pars(18) = DeltaYt;
   pars(19) = DeltaOS;
   pars(20) = Qmatch;
   pars(21) = DeltaLambda3L;
   pars(22) = ThreeLoopAtAsAs;
   pars(23) = msq2(0,0);
   pars(24) = msq2(0,1);
   pars(25) = msq2(0,2);
   pars(26) = msq2(1,0);
   pars(27) = msq2(1,1);
   pars(28) = msq2(1,2);
   pars(29) = msq2(2,0);
   pars(30) = msq2(2,1);
   pars(31) = msq2(2,2);
   pars(32) = msu2(0,0);
   pars(33) = msu2(0,1);
   pars(34) = msu2(0,2);
   pars(35) = msu2(1,0);
   pars(36) = msu2(1,1);
   pars(37) = msu2(1,2);
   pars(38) = msu2(2,0);
   pars(39) = msu2(2,1);
   pars(40) = msu2(2,2);
   pars(41) = msd2(0,0);
   pars(42) = msd2(0,1);
   pars(43) = msd2(0,2);
   pars(44) = msd2(1,0);
   pars(45) = msd2(1,1);
   pars(46) = msd2(1,2);
   pars(47) = msd2(2,0);
   pars(48) = msd2(2,1);
   pars(49) = msd2(2,2);
   pars(50) = msl2(0,0);
   pars(51) = msl2(0,1);
   pars(52) = msl2(0,2);
   pars(53) = msl2(1,0);
   pars(54) = msl2(1,1);
   pars(55) = msl2(1,2);
   pars(56) = msl2(2,0);
   pars(57) = msl2(2,1);
   pars(58) = msl2(2,2);
   pars(59) = mse2(0,0);
   pars(60) = mse2(0,1);
   pars(61) = mse2(0,2);
   pars(62) = mse2(1,0);
   pars(63) = mse2(1,1);
   pars(64) = mse2(1,2);
   pars(65) = mse2(2,0);
   pars(66) = mse2(2,1);
   pars(67) = mse2(2,2);

   return pars;
}

void HSSUSY_input_parameters::set(const Eigen::ArrayXd& pars)
{
   MSUSY = pars(0);
   M1Input = pars(1);
   M2Input = pars(2);
   M3Input = pars(3);
   MuInput = pars(4);
   mAInput = pars(5);
   MEWSB = pars(6);
   AtInput = pars(7);
   AbInput = pars(8);
   AtauInput = pars(9);
   TanBeta = pars(10);
   LambdaLoopOrder = pars(11);
   TwoLoopAtAs = pars(12);
   TwoLoopAbAs = pars(13);
   TwoLoopAtAb = pars(14);
   TwoLoopAtauAtau = pars(15);
   TwoLoopAtAt = pars(16);
   DeltaEFT = pars(17);
   DeltaYt = pars(18);
   DeltaOS = pars(19);
   Qmatch = pars(20);
   DeltaLambda3L = pars(21);
   ThreeLoopAtAsAs = pars(22);
   msq2(0,0) = pars(23);
   msq2(0,1) = pars(24);
   msq2(0,2) = pars(25);
   msq2(1,0) = pars(26);
   msq2(1,1) = pars(27);
   msq2(1,2) = pars(28);
   msq2(2,0) = pars(29);
   msq2(2,1) = pars(30);
   msq2(2,2) = pars(31);
   msu2(0,0) = pars(32);
   msu2(0,1) = pars(33);
   msu2(0,2) = pars(34);
   msu2(1,0) = pars(35);
   msu2(1,1) = pars(36);
   msu2(1,2) = pars(37);
   msu2(2,0) = pars(38);
   msu2(2,1) = pars(39);
   msu2(2,2) = pars(40);
   msd2(0,0) = pars(41);
   msd2(0,1) = pars(42);
   msd2(0,2) = pars(43);
   msd2(1,0) = pars(44);
   msd2(1,1) = pars(45);
   msd2(1,2) = pars(46);
   msd2(2,0) = pars(47);
   msd2(2,1) = pars(48);
   msd2(2,2) = pars(49);
   msl2(0,0) = pars(50);
   msl2(0,1) = pars(51);
   msl2(0,2) = pars(52);
   msl2(1,0) = pars(53);
   msl2(1,1) = pars(54);
   msl2(1,2) = pars(55);
   msl2(2,0) = pars(56);
   msl2(2,1) = pars(57);
   msl2(2,2) = pars(58);
   mse2(0,0) = pars(59);
   mse2(0,1) = pars(60);
   mse2(0,2) = pars(61);
   mse2(1,0) = pars(62);
   mse2(1,1) = pars(63);
   mse2(1,2) = pars(64);
   mse2(2,0) = pars(65);
   mse2(2,1) = pars(66);
   mse2(2,2) = pars(67);

}

std::ostream& operator<<(std::ostream& ostr, const HSSUSY_input_parameters& input)
{
   ostr << "MSUSY = " << INPUT(MSUSY) << ", ";
   ostr << "M1Input = " << INPUT(M1Input) << ", ";
   ostr << "M2Input = " << INPUT(M2Input) << ", ";
   ostr << "M3Input = " << INPUT(M3Input) << ", ";
   ostr << "MuInput = " << INPUT(MuInput) << ", ";
   ostr << "mAInput = " << INPUT(mAInput) << ", ";
   ostr << "MEWSB = " << INPUT(MEWSB) << ", ";
   ostr << "AtInput = " << INPUT(AtInput) << ", ";
   ostr << "AbInput = " << INPUT(AbInput) << ", ";
   ostr << "AtauInput = " << INPUT(AtauInput) << ", ";
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "LambdaLoopOrder = " << INPUT(LambdaLoopOrder) << ", ";
   ostr << "TwoLoopAtAs = " << INPUT(TwoLoopAtAs) << ", ";
   ostr << "TwoLoopAbAs = " << INPUT(TwoLoopAbAs) << ", ";
   ostr << "TwoLoopAtAb = " << INPUT(TwoLoopAtAb) << ", ";
   ostr << "TwoLoopAtauAtau = " << INPUT(TwoLoopAtauAtau) << ", ";
   ostr << "TwoLoopAtAt = " << INPUT(TwoLoopAtAt) << ", ";
   ostr << "DeltaEFT = " << INPUT(DeltaEFT) << ", ";
   ostr << "DeltaYt = " << INPUT(DeltaYt) << ", ";
   ostr << "DeltaOS = " << INPUT(DeltaOS) << ", ";
   ostr << "Qmatch = " << INPUT(Qmatch) << ", ";
   ostr << "DeltaLambda3L = " << INPUT(DeltaLambda3L) << ", ";
   ostr << "ThreeLoopAtAsAs = " << INPUT(ThreeLoopAtAsAs) << ", ";
   ostr << "msq2 = " << INPUT(msq2) << ", ";
   ostr << "msu2 = " << INPUT(msu2) << ", ";
   ostr << "msd2 = " << INPUT(msd2) << ", ";
   ostr << "msl2 = " << INPUT(msl2) << ", ";
   ostr << "mse2 = " << INPUT(mse2) << ", ";

   return ostr;
}

} // namespace flexiblesusy
