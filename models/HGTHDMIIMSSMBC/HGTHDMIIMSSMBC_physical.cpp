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

// File generated at Mon 5 Mar 2018 17:14:19

#include "HGTHDMIIMSSMBC_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

void HGTHDMIIMSSMBC_physical::clear()
{
   MVG = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MGlu = 0.;
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
   MChi = Eigen::Matrix<double,4,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MCha = Eigen::Matrix<double,2,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;

}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void HGTHDMIIMSSMBC_physical::convert_to_hk()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void HGTHDMIIMSSMBC_physical::convert_to_slha()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

Eigen::ArrayXd HGTHDMIIMSSMBC_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(201);

   pars(29) = ZH(0,0);
   pars(30) = ZH(0,1);
   pars(31) = ZH(1,0);
   pars(32) = ZH(1,1);
   pars(33) = ZA(0,0);
   pars(34) = ZA(0,1);
   pars(35) = ZA(1,0);
   pars(36) = ZA(1,1);
   pars(37) = ZP(0,0);
   pars(38) = ZP(0,1);
   pars(39) = ZP(1,0);
   pars(40) = ZP(1,1);
   pars(41) = Re(Vd(0,0));
   pars(42) = Im(Vd(0,0));
   pars(43) = Re(Vd(0,1));
   pars(44) = Im(Vd(0,1));
   pars(45) = Re(Vd(0,2));
   pars(46) = Im(Vd(0,2));
   pars(47) = Re(Vd(1,0));
   pars(48) = Im(Vd(1,0));
   pars(49) = Re(Vd(1,1));
   pars(50) = Im(Vd(1,1));
   pars(51) = Re(Vd(1,2));
   pars(52) = Im(Vd(1,2));
   pars(53) = Re(Vd(2,0));
   pars(54) = Im(Vd(2,0));
   pars(55) = Re(Vd(2,1));
   pars(56) = Im(Vd(2,1));
   pars(57) = Re(Vd(2,2));
   pars(58) = Im(Vd(2,2));
   pars(59) = Re(Ud(0,0));
   pars(60) = Im(Ud(0,0));
   pars(61) = Re(Ud(0,1));
   pars(62) = Im(Ud(0,1));
   pars(63) = Re(Ud(0,2));
   pars(64) = Im(Ud(0,2));
   pars(65) = Re(Ud(1,0));
   pars(66) = Im(Ud(1,0));
   pars(67) = Re(Ud(1,1));
   pars(68) = Im(Ud(1,1));
   pars(69) = Re(Ud(1,2));
   pars(70) = Im(Ud(1,2));
   pars(71) = Re(Ud(2,0));
   pars(72) = Im(Ud(2,0));
   pars(73) = Re(Ud(2,1));
   pars(74) = Im(Ud(2,1));
   pars(75) = Re(Ud(2,2));
   pars(76) = Im(Ud(2,2));
   pars(77) = Re(Vu(0,0));
   pars(78) = Im(Vu(0,0));
   pars(79) = Re(Vu(0,1));
   pars(80) = Im(Vu(0,1));
   pars(81) = Re(Vu(0,2));
   pars(82) = Im(Vu(0,2));
   pars(83) = Re(Vu(1,0));
   pars(84) = Im(Vu(1,0));
   pars(85) = Re(Vu(1,1));
   pars(86) = Im(Vu(1,1));
   pars(87) = Re(Vu(1,2));
   pars(88) = Im(Vu(1,2));
   pars(89) = Re(Vu(2,0));
   pars(90) = Im(Vu(2,0));
   pars(91) = Re(Vu(2,1));
   pars(92) = Im(Vu(2,1));
   pars(93) = Re(Vu(2,2));
   pars(94) = Im(Vu(2,2));
   pars(95) = Re(Uu(0,0));
   pars(96) = Im(Uu(0,0));
   pars(97) = Re(Uu(0,1));
   pars(98) = Im(Uu(0,1));
   pars(99) = Re(Uu(0,2));
   pars(100) = Im(Uu(0,2));
   pars(101) = Re(Uu(1,0));
   pars(102) = Im(Uu(1,0));
   pars(103) = Re(Uu(1,1));
   pars(104) = Im(Uu(1,1));
   pars(105) = Re(Uu(1,2));
   pars(106) = Im(Uu(1,2));
   pars(107) = Re(Uu(2,0));
   pars(108) = Im(Uu(2,0));
   pars(109) = Re(Uu(2,1));
   pars(110) = Im(Uu(2,1));
   pars(111) = Re(Uu(2,2));
   pars(112) = Im(Uu(2,2));
   pars(113) = Re(Ve(0,0));
   pars(114) = Im(Ve(0,0));
   pars(115) = Re(Ve(0,1));
   pars(116) = Im(Ve(0,1));
   pars(117) = Re(Ve(0,2));
   pars(118) = Im(Ve(0,2));
   pars(119) = Re(Ve(1,0));
   pars(120) = Im(Ve(1,0));
   pars(121) = Re(Ve(1,1));
   pars(122) = Im(Ve(1,1));
   pars(123) = Re(Ve(1,2));
   pars(124) = Im(Ve(1,2));
   pars(125) = Re(Ve(2,0));
   pars(126) = Im(Ve(2,0));
   pars(127) = Re(Ve(2,1));
   pars(128) = Im(Ve(2,1));
   pars(129) = Re(Ve(2,2));
   pars(130) = Im(Ve(2,2));
   pars(131) = Re(Ue(0,0));
   pars(132) = Im(Ue(0,0));
   pars(133) = Re(Ue(0,1));
   pars(134) = Im(Ue(0,1));
   pars(135) = Re(Ue(0,2));
   pars(136) = Im(Ue(0,2));
   pars(137) = Re(Ue(1,0));
   pars(138) = Im(Ue(1,0));
   pars(139) = Re(Ue(1,1));
   pars(140) = Im(Ue(1,1));
   pars(141) = Re(Ue(1,2));
   pars(142) = Im(Ue(1,2));
   pars(143) = Re(Ue(2,0));
   pars(144) = Im(Ue(2,0));
   pars(145) = Re(Ue(2,1));
   pars(146) = Im(Ue(2,1));
   pars(147) = Re(Ue(2,2));
   pars(148) = Im(Ue(2,2));
   pars(149) = Re(ZN(0,0));
   pars(150) = Im(ZN(0,0));
   pars(151) = Re(ZN(0,1));
   pars(152) = Im(ZN(0,1));
   pars(153) = Re(ZN(0,2));
   pars(154) = Im(ZN(0,2));
   pars(155) = Re(ZN(0,3));
   pars(156) = Im(ZN(0,3));
   pars(157) = Re(ZN(1,0));
   pars(158) = Im(ZN(1,0));
   pars(159) = Re(ZN(1,1));
   pars(160) = Im(ZN(1,1));
   pars(161) = Re(ZN(1,2));
   pars(162) = Im(ZN(1,2));
   pars(163) = Re(ZN(1,3));
   pars(164) = Im(ZN(1,3));
   pars(165) = Re(ZN(2,0));
   pars(166) = Im(ZN(2,0));
   pars(167) = Re(ZN(2,1));
   pars(168) = Im(ZN(2,1));
   pars(169) = Re(ZN(2,2));
   pars(170) = Im(ZN(2,2));
   pars(171) = Re(ZN(2,3));
   pars(172) = Im(ZN(2,3));
   pars(173) = Re(ZN(3,0));
   pars(174) = Im(ZN(3,0));
   pars(175) = Re(ZN(3,1));
   pars(176) = Im(ZN(3,1));
   pars(177) = Re(ZN(3,2));
   pars(178) = Im(ZN(3,2));
   pars(179) = Re(ZN(3,3));
   pars(180) = Im(ZN(3,3));
   pars(181) = Re(UM(0,0));
   pars(182) = Im(UM(0,0));
   pars(183) = Re(UM(0,1));
   pars(184) = Im(UM(0,1));
   pars(185) = Re(UM(1,0));
   pars(186) = Im(UM(1,0));
   pars(187) = Re(UM(1,1));
   pars(188) = Im(UM(1,1));
   pars(189) = Re(UP(0,0));
   pars(190) = Im(UP(0,0));
   pars(191) = Re(UP(0,1));
   pars(192) = Im(UP(0,1));
   pars(193) = Re(UP(1,0));
   pars(194) = Im(UP(1,0));
   pars(195) = Re(UP(1,1));
   pars(196) = Im(UP(1,1));
   pars(197) = ZZ(0,0);
   pars(198) = ZZ(0,1);
   pars(199) = ZZ(1,0);
   pars(200) = ZZ(1,1);


   return pars;
}

void HGTHDMIIMSSMBC_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

   ZH(0,0) = pars(29);
   ZH(0,1) = pars(30);
   ZH(1,0) = pars(31);
   ZH(1,1) = pars(32);
   ZA(0,0) = pars(33);
   ZA(0,1) = pars(34);
   ZA(1,0) = pars(35);
   ZA(1,1) = pars(36);
   ZP(0,0) = pars(37);
   ZP(0,1) = pars(38);
   ZP(1,0) = pars(39);
   ZP(1,1) = pars(40);
   Vd(0,0) = std::complex<double>(pars(41), pars(42));
   Vd(0,1) = std::complex<double>(pars(43), pars(44));
   Vd(0,2) = std::complex<double>(pars(45), pars(46));
   Vd(1,0) = std::complex<double>(pars(47), pars(48));
   Vd(1,1) = std::complex<double>(pars(49), pars(50));
   Vd(1,2) = std::complex<double>(pars(51), pars(52));
   Vd(2,0) = std::complex<double>(pars(53), pars(54));
   Vd(2,1) = std::complex<double>(pars(55), pars(56));
   Vd(2,2) = std::complex<double>(pars(57), pars(58));
   Ud(0,0) = std::complex<double>(pars(59), pars(60));
   Ud(0,1) = std::complex<double>(pars(61), pars(62));
   Ud(0,2) = std::complex<double>(pars(63), pars(64));
   Ud(1,0) = std::complex<double>(pars(65), pars(66));
   Ud(1,1) = std::complex<double>(pars(67), pars(68));
   Ud(1,2) = std::complex<double>(pars(69), pars(70));
   Ud(2,0) = std::complex<double>(pars(71), pars(72));
   Ud(2,1) = std::complex<double>(pars(73), pars(74));
   Ud(2,2) = std::complex<double>(pars(75), pars(76));
   Vu(0,0) = std::complex<double>(pars(77), pars(78));
   Vu(0,1) = std::complex<double>(pars(79), pars(80));
   Vu(0,2) = std::complex<double>(pars(81), pars(82));
   Vu(1,0) = std::complex<double>(pars(83), pars(84));
   Vu(1,1) = std::complex<double>(pars(85), pars(86));
   Vu(1,2) = std::complex<double>(pars(87), pars(88));
   Vu(2,0) = std::complex<double>(pars(89), pars(90));
   Vu(2,1) = std::complex<double>(pars(91), pars(92));
   Vu(2,2) = std::complex<double>(pars(93), pars(94));
   Uu(0,0) = std::complex<double>(pars(95), pars(96));
   Uu(0,1) = std::complex<double>(pars(97), pars(98));
   Uu(0,2) = std::complex<double>(pars(99), pars(100));
   Uu(1,0) = std::complex<double>(pars(101), pars(102));
   Uu(1,1) = std::complex<double>(pars(103), pars(104));
   Uu(1,2) = std::complex<double>(pars(105), pars(106));
   Uu(2,0) = std::complex<double>(pars(107), pars(108));
   Uu(2,1) = std::complex<double>(pars(109), pars(110));
   Uu(2,2) = std::complex<double>(pars(111), pars(112));
   Ve(0,0) = std::complex<double>(pars(113), pars(114));
   Ve(0,1) = std::complex<double>(pars(115), pars(116));
   Ve(0,2) = std::complex<double>(pars(117), pars(118));
   Ve(1,0) = std::complex<double>(pars(119), pars(120));
   Ve(1,1) = std::complex<double>(pars(121), pars(122));
   Ve(1,2) = std::complex<double>(pars(123), pars(124));
   Ve(2,0) = std::complex<double>(pars(125), pars(126));
   Ve(2,1) = std::complex<double>(pars(127), pars(128));
   Ve(2,2) = std::complex<double>(pars(129), pars(130));
   Ue(0,0) = std::complex<double>(pars(131), pars(132));
   Ue(0,1) = std::complex<double>(pars(133), pars(134));
   Ue(0,2) = std::complex<double>(pars(135), pars(136));
   Ue(1,0) = std::complex<double>(pars(137), pars(138));
   Ue(1,1) = std::complex<double>(pars(139), pars(140));
   Ue(1,2) = std::complex<double>(pars(141), pars(142));
   Ue(2,0) = std::complex<double>(pars(143), pars(144));
   Ue(2,1) = std::complex<double>(pars(145), pars(146));
   Ue(2,2) = std::complex<double>(pars(147), pars(148));
   ZN(0,0) = std::complex<double>(pars(149), pars(150));
   ZN(0,1) = std::complex<double>(pars(151), pars(152));
   ZN(0,2) = std::complex<double>(pars(153), pars(154));
   ZN(0,3) = std::complex<double>(pars(155), pars(156));
   ZN(1,0) = std::complex<double>(pars(157), pars(158));
   ZN(1,1) = std::complex<double>(pars(159), pars(160));
   ZN(1,2) = std::complex<double>(pars(161), pars(162));
   ZN(1,3) = std::complex<double>(pars(163), pars(164));
   ZN(2,0) = std::complex<double>(pars(165), pars(166));
   ZN(2,1) = std::complex<double>(pars(167), pars(168));
   ZN(2,2) = std::complex<double>(pars(169), pars(170));
   ZN(2,3) = std::complex<double>(pars(171), pars(172));
   ZN(3,0) = std::complex<double>(pars(173), pars(174));
   ZN(3,1) = std::complex<double>(pars(175), pars(176));
   ZN(3,2) = std::complex<double>(pars(177), pars(178));
   ZN(3,3) = std::complex<double>(pars(179), pars(180));
   UM(0,0) = std::complex<double>(pars(181), pars(182));
   UM(0,1) = std::complex<double>(pars(183), pars(184));
   UM(1,0) = std::complex<double>(pars(185), pars(186));
   UM(1,1) = std::complex<double>(pars(187), pars(188));
   UP(0,0) = std::complex<double>(pars(189), pars(190));
   UP(0,1) = std::complex<double>(pars(191), pars(192));
   UP(1,0) = std::complex<double>(pars(193), pars(194));
   UP(1,1) = std::complex<double>(pars(195), pars(196));
   ZZ(0,0) = pars(197);
   ZZ(0,1) = pars(198);
   ZZ(1,0) = pars(199);
   ZZ(1,1) = pars(200);

}

Eigen::ArrayXd HGTHDMIIMSSMBC_physical::get_masses() const
{
   Eigen::ArrayXd pars(29);

   pars(0) = MVG;
   pars(1) = MFv(0);
   pars(2) = MFv(1);
   pars(3) = MFv(2);
   pars(4) = MGlu;
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
   pars(22) = MChi(2);
   pars(23) = MChi(3);
   pars(24) = MCha(0);
   pars(25) = MCha(1);
   pars(26) = MVWm;
   pars(27) = MVP;
   pars(28) = MVZ;

   return pars;
}

void HGTHDMIIMSSMBC_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MFv(0) = pars(1);
   MFv(1) = pars(2);
   MFv(2) = pars(3);
   MGlu = pars(4);
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
   MChi(2) = pars(22);
   MChi(3) = pars(23);
   MCha(0) = pars(24);
   MCha(1) = pars(25);
   MVWm = pars(26);
   MVP = pars(27);
   MVZ = pars(28);

}

void HGTHDMIIMSSMBC_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHm = " << MHm.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
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
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
   ostr << "ZZ = " << ZZ << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const HGTHDMIIMSSMBC_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
