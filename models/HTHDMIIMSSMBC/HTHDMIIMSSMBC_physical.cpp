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

// File generated at Sun 4 Aug 2019 18:59:26

#include "HTHDMIIMSSMBC_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

void HTHDMIIMSSMBC_physical::clear()
{
   MVG = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MCha = 0.;
   Mhh = Eigen::Matrix<double,2,1>::Zero();
   ZH = Eigen::Matrix<double,2,2>::Zero();
   MAh = Eigen::Matrix<double,2,1>::Zero();
   ZA = Eigen::Matrix<double,2,2>::Zero();
   MHm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MChi = Eigen::Matrix<double,2,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;

}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void HTHDMIIMSSMBC_physical::convert_to_hk()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void HTHDMIIMSSMBC_physical::convert_to_slha()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

Eigen::ArrayXd HTHDMIIMSSMBC_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(157);

   pars(25) = ZH(0,0);
   pars(26) = ZH(0,1);
   pars(27) = ZH(1,0);
   pars(28) = ZH(1,1);
   pars(29) = ZA(0,0);
   pars(30) = ZA(0,1);
   pars(31) = ZA(1,0);
   pars(32) = ZA(1,1);
   pars(33) = ZP(0,0);
   pars(34) = ZP(0,1);
   pars(35) = ZP(1,0);
   pars(36) = ZP(1,1);
   pars(37) = Re(Vd(0,0));
   pars(38) = Im(Vd(0,0));
   pars(39) = Re(Vd(0,1));
   pars(40) = Im(Vd(0,1));
   pars(41) = Re(Vd(0,2));
   pars(42) = Im(Vd(0,2));
   pars(43) = Re(Vd(1,0));
   pars(44) = Im(Vd(1,0));
   pars(45) = Re(Vd(1,1));
   pars(46) = Im(Vd(1,1));
   pars(47) = Re(Vd(1,2));
   pars(48) = Im(Vd(1,2));
   pars(49) = Re(Vd(2,0));
   pars(50) = Im(Vd(2,0));
   pars(51) = Re(Vd(2,1));
   pars(52) = Im(Vd(2,1));
   pars(53) = Re(Vd(2,2));
   pars(54) = Im(Vd(2,2));
   pars(55) = Re(Ud(0,0));
   pars(56) = Im(Ud(0,0));
   pars(57) = Re(Ud(0,1));
   pars(58) = Im(Ud(0,1));
   pars(59) = Re(Ud(0,2));
   pars(60) = Im(Ud(0,2));
   pars(61) = Re(Ud(1,0));
   pars(62) = Im(Ud(1,0));
   pars(63) = Re(Ud(1,1));
   pars(64) = Im(Ud(1,1));
   pars(65) = Re(Ud(1,2));
   pars(66) = Im(Ud(1,2));
   pars(67) = Re(Ud(2,0));
   pars(68) = Im(Ud(2,0));
   pars(69) = Re(Ud(2,1));
   pars(70) = Im(Ud(2,1));
   pars(71) = Re(Ud(2,2));
   pars(72) = Im(Ud(2,2));
   pars(73) = Re(Vu(0,0));
   pars(74) = Im(Vu(0,0));
   pars(75) = Re(Vu(0,1));
   pars(76) = Im(Vu(0,1));
   pars(77) = Re(Vu(0,2));
   pars(78) = Im(Vu(0,2));
   pars(79) = Re(Vu(1,0));
   pars(80) = Im(Vu(1,0));
   pars(81) = Re(Vu(1,1));
   pars(82) = Im(Vu(1,1));
   pars(83) = Re(Vu(1,2));
   pars(84) = Im(Vu(1,2));
   pars(85) = Re(Vu(2,0));
   pars(86) = Im(Vu(2,0));
   pars(87) = Re(Vu(2,1));
   pars(88) = Im(Vu(2,1));
   pars(89) = Re(Vu(2,2));
   pars(90) = Im(Vu(2,2));
   pars(91) = Re(Uu(0,0));
   pars(92) = Im(Uu(0,0));
   pars(93) = Re(Uu(0,1));
   pars(94) = Im(Uu(0,1));
   pars(95) = Re(Uu(0,2));
   pars(96) = Im(Uu(0,2));
   pars(97) = Re(Uu(1,0));
   pars(98) = Im(Uu(1,0));
   pars(99) = Re(Uu(1,1));
   pars(100) = Im(Uu(1,1));
   pars(101) = Re(Uu(1,2));
   pars(102) = Im(Uu(1,2));
   pars(103) = Re(Uu(2,0));
   pars(104) = Im(Uu(2,0));
   pars(105) = Re(Uu(2,1));
   pars(106) = Im(Uu(2,1));
   pars(107) = Re(Uu(2,2));
   pars(108) = Im(Uu(2,2));
   pars(109) = Re(Ve(0,0));
   pars(110) = Im(Ve(0,0));
   pars(111) = Re(Ve(0,1));
   pars(112) = Im(Ve(0,1));
   pars(113) = Re(Ve(0,2));
   pars(114) = Im(Ve(0,2));
   pars(115) = Re(Ve(1,0));
   pars(116) = Im(Ve(1,0));
   pars(117) = Re(Ve(1,1));
   pars(118) = Im(Ve(1,1));
   pars(119) = Re(Ve(1,2));
   pars(120) = Im(Ve(1,2));
   pars(121) = Re(Ve(2,0));
   pars(122) = Im(Ve(2,0));
   pars(123) = Re(Ve(2,1));
   pars(124) = Im(Ve(2,1));
   pars(125) = Re(Ve(2,2));
   pars(126) = Im(Ve(2,2));
   pars(127) = Re(Ue(0,0));
   pars(128) = Im(Ue(0,0));
   pars(129) = Re(Ue(0,1));
   pars(130) = Im(Ue(0,1));
   pars(131) = Re(Ue(0,2));
   pars(132) = Im(Ue(0,2));
   pars(133) = Re(Ue(1,0));
   pars(134) = Im(Ue(1,0));
   pars(135) = Re(Ue(1,1));
   pars(136) = Im(Ue(1,1));
   pars(137) = Re(Ue(1,2));
   pars(138) = Im(Ue(1,2));
   pars(139) = Re(Ue(2,0));
   pars(140) = Im(Ue(2,0));
   pars(141) = Re(Ue(2,1));
   pars(142) = Im(Ue(2,1));
   pars(143) = Re(Ue(2,2));
   pars(144) = Im(Ue(2,2));
   pars(145) = Re(ZN(0,0));
   pars(146) = Im(ZN(0,0));
   pars(147) = Re(ZN(0,1));
   pars(148) = Im(ZN(0,1));
   pars(149) = Re(ZN(1,0));
   pars(150) = Im(ZN(1,0));
   pars(151) = Re(ZN(1,1));
   pars(152) = Im(ZN(1,1));
   pars(153) = ZZ(0,0);
   pars(154) = ZZ(0,1);
   pars(155) = ZZ(1,0);
   pars(156) = ZZ(1,1);


   return pars;
}

void HTHDMIIMSSMBC_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

   ZH(0,0) = pars(25);
   ZH(0,1) = pars(26);
   ZH(1,0) = pars(27);
   ZH(1,1) = pars(28);
   ZA(0,0) = pars(29);
   ZA(0,1) = pars(30);
   ZA(1,0) = pars(31);
   ZA(1,1) = pars(32);
   ZP(0,0) = pars(33);
   ZP(0,1) = pars(34);
   ZP(1,0) = pars(35);
   ZP(1,1) = pars(36);
   Vd(0,0) = std::complex<double>(pars(37), pars(38));
   Vd(0,1) = std::complex<double>(pars(39), pars(40));
   Vd(0,2) = std::complex<double>(pars(41), pars(42));
   Vd(1,0) = std::complex<double>(pars(43), pars(44));
   Vd(1,1) = std::complex<double>(pars(45), pars(46));
   Vd(1,2) = std::complex<double>(pars(47), pars(48));
   Vd(2,0) = std::complex<double>(pars(49), pars(50));
   Vd(2,1) = std::complex<double>(pars(51), pars(52));
   Vd(2,2) = std::complex<double>(pars(53), pars(54));
   Ud(0,0) = std::complex<double>(pars(55), pars(56));
   Ud(0,1) = std::complex<double>(pars(57), pars(58));
   Ud(0,2) = std::complex<double>(pars(59), pars(60));
   Ud(1,0) = std::complex<double>(pars(61), pars(62));
   Ud(1,1) = std::complex<double>(pars(63), pars(64));
   Ud(1,2) = std::complex<double>(pars(65), pars(66));
   Ud(2,0) = std::complex<double>(pars(67), pars(68));
   Ud(2,1) = std::complex<double>(pars(69), pars(70));
   Ud(2,2) = std::complex<double>(pars(71), pars(72));
   Vu(0,0) = std::complex<double>(pars(73), pars(74));
   Vu(0,1) = std::complex<double>(pars(75), pars(76));
   Vu(0,2) = std::complex<double>(pars(77), pars(78));
   Vu(1,0) = std::complex<double>(pars(79), pars(80));
   Vu(1,1) = std::complex<double>(pars(81), pars(82));
   Vu(1,2) = std::complex<double>(pars(83), pars(84));
   Vu(2,0) = std::complex<double>(pars(85), pars(86));
   Vu(2,1) = std::complex<double>(pars(87), pars(88));
   Vu(2,2) = std::complex<double>(pars(89), pars(90));
   Uu(0,0) = std::complex<double>(pars(91), pars(92));
   Uu(0,1) = std::complex<double>(pars(93), pars(94));
   Uu(0,2) = std::complex<double>(pars(95), pars(96));
   Uu(1,0) = std::complex<double>(pars(97), pars(98));
   Uu(1,1) = std::complex<double>(pars(99), pars(100));
   Uu(1,2) = std::complex<double>(pars(101), pars(102));
   Uu(2,0) = std::complex<double>(pars(103), pars(104));
   Uu(2,1) = std::complex<double>(pars(105), pars(106));
   Uu(2,2) = std::complex<double>(pars(107), pars(108));
   Ve(0,0) = std::complex<double>(pars(109), pars(110));
   Ve(0,1) = std::complex<double>(pars(111), pars(112));
   Ve(0,2) = std::complex<double>(pars(113), pars(114));
   Ve(1,0) = std::complex<double>(pars(115), pars(116));
   Ve(1,1) = std::complex<double>(pars(117), pars(118));
   Ve(1,2) = std::complex<double>(pars(119), pars(120));
   Ve(2,0) = std::complex<double>(pars(121), pars(122));
   Ve(2,1) = std::complex<double>(pars(123), pars(124));
   Ve(2,2) = std::complex<double>(pars(125), pars(126));
   Ue(0,0) = std::complex<double>(pars(127), pars(128));
   Ue(0,1) = std::complex<double>(pars(129), pars(130));
   Ue(0,2) = std::complex<double>(pars(131), pars(132));
   Ue(1,0) = std::complex<double>(pars(133), pars(134));
   Ue(1,1) = std::complex<double>(pars(135), pars(136));
   Ue(1,2) = std::complex<double>(pars(137), pars(138));
   Ue(2,0) = std::complex<double>(pars(139), pars(140));
   Ue(2,1) = std::complex<double>(pars(141), pars(142));
   Ue(2,2) = std::complex<double>(pars(143), pars(144));
   ZN(0,0) = std::complex<double>(pars(145), pars(146));
   ZN(0,1) = std::complex<double>(pars(147), pars(148));
   ZN(1,0) = std::complex<double>(pars(149), pars(150));
   ZN(1,1) = std::complex<double>(pars(151), pars(152));
   ZZ(0,0) = pars(153);
   ZZ(0,1) = pars(154);
   ZZ(1,0) = pars(155);
   ZZ(1,1) = pars(156);

}

Eigen::ArrayXd HTHDMIIMSSMBC_physical::get_masses() const
{
   Eigen::ArrayXd pars(25);

   pars(0) = MVG;
   pars(1) = MFv(0);
   pars(2) = MFv(1);
   pars(3) = MFv(2);
   pars(4) = MCha;
   pars(5) = Mhh(0);
   pars(6) = Mhh(1);
   pars(7) = MAh(0);
   pars(8) = MAh(1);
   pars(9) = MHm(0);
   pars(10) = MHm(1);
   pars(11) = MFd(0);
   pars(12) = MFd(1);
   pars(13) = MFd(2);
   pars(14) = MFu(0);
   pars(15) = MFu(1);
   pars(16) = MFu(2);
   pars(17) = MFe(0);
   pars(18) = MFe(1);
   pars(19) = MFe(2);
   pars(20) = MChi(0);
   pars(21) = MChi(1);
   pars(22) = MVWm;
   pars(23) = MVP;
   pars(24) = MVZ;

   return pars;
}

void HTHDMIIMSSMBC_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MFv(0) = pars(1);
   MFv(1) = pars(2);
   MFv(2) = pars(3);
   MCha = pars(4);
   Mhh(0) = pars(5);
   Mhh(1) = pars(6);
   MAh(0) = pars(7);
   MAh(1) = pars(8);
   MHm(0) = pars(9);
   MHm(1) = pars(10);
   MFd(0) = pars(11);
   MFd(1) = pars(12);
   MFd(2) = pars(13);
   MFu(0) = pars(14);
   MFu(1) = pars(15);
   MFu(2) = pars(16);
   MFe(0) = pars(17);
   MFe(1) = pars(18);
   MFe(2) = pars(19);
   MChi(0) = pars(20);
   MChi(1) = pars(21);
   MVWm = pars(22);
   MVP = pars(23);
   MVZ = pars(24);

}

void HTHDMIIMSSMBC_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MCha = " << MCha << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHm = " << MHm.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';
   ostr << "MVP = " << MVP << '\n';
   ostr << "MVZ = " << MVZ << '\n';

   ostr << "----------------------------------------\n"
           "pole mass mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "ZZ = " << ZZ << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const HTHDMIIMSSMBC_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
