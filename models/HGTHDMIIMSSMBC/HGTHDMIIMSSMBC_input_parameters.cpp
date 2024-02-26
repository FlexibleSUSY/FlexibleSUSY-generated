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


#include "HGTHDMIIMSSMBC_input_parameters.hpp"
#include "wrappers.hpp"

#define INPUT(p) input.p

namespace flexiblesusy {

Eigen::ArrayXd HGTHDMIIMSSMBC_input_parameters::get() const
{
   Eigen::ArrayXd pars(51);

   pars(0) = TanBeta;
   pars(1) = MSUSY;
   pars(2) = MEWSB;
   pars(3) = MuInput;
   pars(4) = M1Input;
   pars(5) = M2Input;
   pars(6) = M3Input;
   pars(7) = MAInput;
   pars(8) = LambdaLoopOrder;
   pars(9) = AeInput(0,0);
   pars(10) = AeInput(0,1);
   pars(11) = AeInput(0,2);
   pars(12) = AeInput(1,0);
   pars(13) = AeInput(1,1);
   pars(14) = AeInput(1,2);
   pars(15) = AeInput(2,0);
   pars(16) = AeInput(2,1);
   pars(17) = AeInput(2,2);
   pars(18) = AdInput(0,0);
   pars(19) = AdInput(0,1);
   pars(20) = AdInput(0,2);
   pars(21) = AdInput(1,0);
   pars(22) = AdInput(1,1);
   pars(23) = AdInput(1,2);
   pars(24) = AdInput(2,0);
   pars(25) = AdInput(2,1);
   pars(26) = AdInput(2,2);
   pars(27) = AuInput(0,0);
   pars(28) = AuInput(0,1);
   pars(29) = AuInput(0,2);
   pars(30) = AuInput(1,0);
   pars(31) = AuInput(1,1);
   pars(32) = AuInput(1,2);
   pars(33) = AuInput(2,0);
   pars(34) = AuInput(2,1);
   pars(35) = AuInput(2,2);
   pars(36) = mslInput(0);
   pars(37) = mslInput(1);
   pars(38) = mslInput(2);
   pars(39) = mseInput(0);
   pars(40) = mseInput(1);
   pars(41) = mseInput(2);
   pars(42) = msqInput(0);
   pars(43) = msqInput(1);
   pars(44) = msqInput(2);
   pars(45) = msdInput(0);
   pars(46) = msdInput(1);
   pars(47) = msdInput(2);
   pars(48) = msuInput(0);
   pars(49) = msuInput(1);
   pars(50) = msuInput(2);

   return pars;
}

void HGTHDMIIMSSMBC_input_parameters::set(const Eigen::ArrayXd& pars)
{
   TanBeta = pars(0);
   MSUSY = pars(1);
   MEWSB = pars(2);
   MuInput = pars(3);
   M1Input = pars(4);
   M2Input = pars(5);
   M3Input = pars(6);
   MAInput = pars(7);
   LambdaLoopOrder = pars(8);
   AeInput(0,0) = pars(9);
   AeInput(0,1) = pars(10);
   AeInput(0,2) = pars(11);
   AeInput(1,0) = pars(12);
   AeInput(1,1) = pars(13);
   AeInput(1,2) = pars(14);
   AeInput(2,0) = pars(15);
   AeInput(2,1) = pars(16);
   AeInput(2,2) = pars(17);
   AdInput(0,0) = pars(18);
   AdInput(0,1) = pars(19);
   AdInput(0,2) = pars(20);
   AdInput(1,0) = pars(21);
   AdInput(1,1) = pars(22);
   AdInput(1,2) = pars(23);
   AdInput(2,0) = pars(24);
   AdInput(2,1) = pars(25);
   AdInput(2,2) = pars(26);
   AuInput(0,0) = pars(27);
   AuInput(0,1) = pars(28);
   AuInput(0,2) = pars(29);
   AuInput(1,0) = pars(30);
   AuInput(1,1) = pars(31);
   AuInput(1,2) = pars(32);
   AuInput(2,0) = pars(33);
   AuInput(2,1) = pars(34);
   AuInput(2,2) = pars(35);
   mslInput(0) = pars(36);
   mslInput(1) = pars(37);
   mslInput(2) = pars(38);
   mseInput(0) = pars(39);
   mseInput(1) = pars(40);
   mseInput(2) = pars(41);
   msqInput(0) = pars(42);
   msqInput(1) = pars(43);
   msqInput(2) = pars(44);
   msdInput(0) = pars(45);
   msdInput(1) = pars(46);
   msdInput(2) = pars(47);
   msuInput(0) = pars(48);
   msuInput(1) = pars(49);
   msuInput(2) = pars(50);

}

std::ostream& operator<<(std::ostream& ostr, const HGTHDMIIMSSMBC_input_parameters& input)
{
   ostr << "TanBeta = " << INPUT(TanBeta) << ", ";
   ostr << "MSUSY = " << INPUT(MSUSY) << ", ";
   ostr << "MEWSB = " << INPUT(MEWSB) << ", ";
   ostr << "MuInput = " << INPUT(MuInput) << ", ";
   ostr << "M1Input = " << INPUT(M1Input) << ", ";
   ostr << "M2Input = " << INPUT(M2Input) << ", ";
   ostr << "M3Input = " << INPUT(M3Input) << ", ";
   ostr << "MAInput = " << INPUT(MAInput) << ", ";
   ostr << "LambdaLoopOrder = " << INPUT(LambdaLoopOrder) << ", ";
   ostr << "AeInput = " << INPUT(AeInput) << ", ";
   ostr << "AdInput = " << INPUT(AdInput) << ", ";
   ostr << "AuInput = " << INPUT(AuInput) << ", ";
   ostr << "mslInput = " << INPUT(mslInput.transpose()) << ", ";
   ostr << "mseInput = " << INPUT(mseInput.transpose()) << ", ";
   ostr << "msqInput = " << INPUT(msqInput.transpose()) << ", ";
   ostr << "msdInput = " << INPUT(msdInput.transpose()) << ", ";
   ostr << "msuInput = " << INPUT(msuInput.transpose()) << ", ";

   return ostr;
}

} // namespace flexiblesusy
