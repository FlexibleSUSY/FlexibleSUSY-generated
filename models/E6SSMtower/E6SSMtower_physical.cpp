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

// File generated at Thu 15 Dec 2016 12:50:22

#include "E6SSMtower_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

E6SSMtower_physical::E6SSMtower_physical()
   :
    MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MChaP(0), MSd(
       Eigen::Array<double,6,1>::Zero()), MSv(Eigen::Array<double,3,1>::Zero()),
       MSu(Eigen::Array<double,6,1>::Zero()), MSe(Eigen::Array<double,6,1>::Zero()
       ), MSDX(Eigen::Array<double,6,1>::Zero()), Mhh(Eigen::Array<double,3,1>
       ::Zero()), MAh(Eigen::Array<double,3,1>::Zero()), MHpm(Eigen::Array<double,
       2,1>::Zero()), MChi(Eigen::Array<double,6,1>::Zero()), MCha(Eigen::Array<
       double,2,1>::Zero()), MFe(Eigen::Array<double,3,1>::Zero()), MFd(
       Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double,3,1>::Zero()),
       MFDX(Eigen::Array<double,3,1>::Zero()), MSHI0(Eigen::Array<double,4,1>
       ::Zero()), MSHIp(Eigen::Array<double,4,1>::Zero()), MChaI(Eigen::Array<
       double,2,1>::Zero()), MChiI(Eigen::Array<double,4,1>::Zero()), MSSI0(
       Eigen::Array<double,2,1>::Zero()), MFSI(Eigen::Array<double,2,1>::Zero()),
       MSHp0(Eigen::Array<double,2,1>::Zero()), MSHpp(Eigen::Array<double,2,1>
       ::Zero()), MChiP(Eigen::Array<double,2,1>::Zero()), MVWm(0), MVP(0), MVZ(0)
       , MVZp(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZDX(Eigen::Matrix<double,6,6>::Zero()), ZH(Eigen::Matrix<double,3
      ,3>::Zero()), ZA(Eigen::Matrix<double,3,3>::Zero()), ZP(Eigen::Matrix<double
      ,2,2>::Zero()), ZN(Eigen::Matrix<std::complex<double>,6,6>::Zero()), UM(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP(Eigen::Matrix<
      std::complex<double>,2,2>::Zero()), ZEL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZER(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZUL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZUR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDXL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDXR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), UHI0(Eigen::Matrix<double,4,4>::Zero()),
      UHIp(Eigen::Matrix<double,4,4>::Zero()), ZMI(Eigen::Matrix<std::complex<
      double>,2,2>::Zero()), ZPI(Eigen::Matrix<std::complex<double>,2,2>::Zero()),
      ZNI(Eigen::Matrix<std::complex<double>,4,4>::Zero()), ZSSI(Eigen::Matrix<
      double,2,2>::Zero()), ZFSI(Eigen::Matrix<std::complex<double>,2,2>::Zero()),
      UHp0(Eigen::Matrix<double,2,2>::Zero()), UHpp(Eigen::Matrix<double,2,2>
      ::Zero()), ZNp(Eigen::Matrix<std::complex<double>,2,2>::Zero()), ZZ(
      Eigen::Matrix<double,3,3>::Zero())

{
}

void E6SSMtower_physical::clear()
{
   MVG = 0.;
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MChaP = 0.;
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   MSDX = Eigen::Matrix<double,6,1>::Zero();
   ZDX = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,6,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,6,6>::Zero();
   MCha = Eigen::Matrix<double,2,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   ZEL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZER = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   ZDL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   ZUL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZUR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFDX = Eigen::Matrix<double,3,1>::Zero();
   ZDXL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDXR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MSHI0 = Eigen::Matrix<double,4,1>::Zero();
   UHI0 = Eigen::Matrix<double,4,4>::Zero();
   MSHIp = Eigen::Matrix<double,4,1>::Zero();
   UHIp = Eigen::Matrix<double,4,4>::Zero();
   MChaI = Eigen::Matrix<double,2,1>::Zero();
   ZMI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   ZPI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MChiI = Eigen::Matrix<double,4,1>::Zero();
   ZNI = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MSSI0 = Eigen::Matrix<double,2,1>::Zero();
   ZSSI = Eigen::Matrix<double,2,2>::Zero();
   MFSI = Eigen::Matrix<double,2,1>::Zero();
   ZFSI = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MSHp0 = Eigen::Matrix<double,2,1>::Zero();
   UHp0 = Eigen::Matrix<double,2,2>::Zero();
   MSHpp = Eigen::Matrix<double,2,1>::Zero();
   UHpp = Eigen::Matrix<double,2,2>::Zero();
   MChiP = Eigen::Matrix<double,2,1>::Zero();
   ZNp = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;
   MVZp = 0.;

}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void E6SSMtower_physical::convert_to_hk()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChiI), LOCALPHYSICAL(ZNI));
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MFSI), LOCALPHYSICAL(ZFSI));
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChiP), LOCALPHYSICAL(ZNp));

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void E6SSMtower_physical::convert_to_slha()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChiI), LOCALPHYSICAL(ZNI));
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MFSI), LOCALPHYSICAL(ZFSI));
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChiP), LOCALPHYSICAL(ZNp));

}

Eigen::ArrayXd E6SSMtower_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(613);

   pars(89) = ZD(0,0);
   pars(90) = ZD(0,1);
   pars(91) = ZD(0,2);
   pars(92) = ZD(0,3);
   pars(93) = ZD(0,4);
   pars(94) = ZD(0,5);
   pars(95) = ZD(1,0);
   pars(96) = ZD(1,1);
   pars(97) = ZD(1,2);
   pars(98) = ZD(1,3);
   pars(99) = ZD(1,4);
   pars(100) = ZD(1,5);
   pars(101) = ZD(2,0);
   pars(102) = ZD(2,1);
   pars(103) = ZD(2,2);
   pars(104) = ZD(2,3);
   pars(105) = ZD(2,4);
   pars(106) = ZD(2,5);
   pars(107) = ZD(3,0);
   pars(108) = ZD(3,1);
   pars(109) = ZD(3,2);
   pars(110) = ZD(3,3);
   pars(111) = ZD(3,4);
   pars(112) = ZD(3,5);
   pars(113) = ZD(4,0);
   pars(114) = ZD(4,1);
   pars(115) = ZD(4,2);
   pars(116) = ZD(4,3);
   pars(117) = ZD(4,4);
   pars(118) = ZD(4,5);
   pars(119) = ZD(5,0);
   pars(120) = ZD(5,1);
   pars(121) = ZD(5,2);
   pars(122) = ZD(5,3);
   pars(123) = ZD(5,4);
   pars(124) = ZD(5,5);
   pars(125) = ZV(0,0);
   pars(126) = ZV(0,1);
   pars(127) = ZV(0,2);
   pars(128) = ZV(1,0);
   pars(129) = ZV(1,1);
   pars(130) = ZV(1,2);
   pars(131) = ZV(2,0);
   pars(132) = ZV(2,1);
   pars(133) = ZV(2,2);
   pars(134) = ZU(0,0);
   pars(135) = ZU(0,1);
   pars(136) = ZU(0,2);
   pars(137) = ZU(0,3);
   pars(138) = ZU(0,4);
   pars(139) = ZU(0,5);
   pars(140) = ZU(1,0);
   pars(141) = ZU(1,1);
   pars(142) = ZU(1,2);
   pars(143) = ZU(1,3);
   pars(144) = ZU(1,4);
   pars(145) = ZU(1,5);
   pars(146) = ZU(2,0);
   pars(147) = ZU(2,1);
   pars(148) = ZU(2,2);
   pars(149) = ZU(2,3);
   pars(150) = ZU(2,4);
   pars(151) = ZU(2,5);
   pars(152) = ZU(3,0);
   pars(153) = ZU(3,1);
   pars(154) = ZU(3,2);
   pars(155) = ZU(3,3);
   pars(156) = ZU(3,4);
   pars(157) = ZU(3,5);
   pars(158) = ZU(4,0);
   pars(159) = ZU(4,1);
   pars(160) = ZU(4,2);
   pars(161) = ZU(4,3);
   pars(162) = ZU(4,4);
   pars(163) = ZU(4,5);
   pars(164) = ZU(5,0);
   pars(165) = ZU(5,1);
   pars(166) = ZU(5,2);
   pars(167) = ZU(5,3);
   pars(168) = ZU(5,4);
   pars(169) = ZU(5,5);
   pars(170) = ZE(0,0);
   pars(171) = ZE(0,1);
   pars(172) = ZE(0,2);
   pars(173) = ZE(0,3);
   pars(174) = ZE(0,4);
   pars(175) = ZE(0,5);
   pars(176) = ZE(1,0);
   pars(177) = ZE(1,1);
   pars(178) = ZE(1,2);
   pars(179) = ZE(1,3);
   pars(180) = ZE(1,4);
   pars(181) = ZE(1,5);
   pars(182) = ZE(2,0);
   pars(183) = ZE(2,1);
   pars(184) = ZE(2,2);
   pars(185) = ZE(2,3);
   pars(186) = ZE(2,4);
   pars(187) = ZE(2,5);
   pars(188) = ZE(3,0);
   pars(189) = ZE(3,1);
   pars(190) = ZE(3,2);
   pars(191) = ZE(3,3);
   pars(192) = ZE(3,4);
   pars(193) = ZE(3,5);
   pars(194) = ZE(4,0);
   pars(195) = ZE(4,1);
   pars(196) = ZE(4,2);
   pars(197) = ZE(4,3);
   pars(198) = ZE(4,4);
   pars(199) = ZE(4,5);
   pars(200) = ZE(5,0);
   pars(201) = ZE(5,1);
   pars(202) = ZE(5,2);
   pars(203) = ZE(5,3);
   pars(204) = ZE(5,4);
   pars(205) = ZE(5,5);
   pars(206) = ZDX(0,0);
   pars(207) = ZDX(0,1);
   pars(208) = ZDX(0,2);
   pars(209) = ZDX(0,3);
   pars(210) = ZDX(0,4);
   pars(211) = ZDX(0,5);
   pars(212) = ZDX(1,0);
   pars(213) = ZDX(1,1);
   pars(214) = ZDX(1,2);
   pars(215) = ZDX(1,3);
   pars(216) = ZDX(1,4);
   pars(217) = ZDX(1,5);
   pars(218) = ZDX(2,0);
   pars(219) = ZDX(2,1);
   pars(220) = ZDX(2,2);
   pars(221) = ZDX(2,3);
   pars(222) = ZDX(2,4);
   pars(223) = ZDX(2,5);
   pars(224) = ZDX(3,0);
   pars(225) = ZDX(3,1);
   pars(226) = ZDX(3,2);
   pars(227) = ZDX(3,3);
   pars(228) = ZDX(3,4);
   pars(229) = ZDX(3,5);
   pars(230) = ZDX(4,0);
   pars(231) = ZDX(4,1);
   pars(232) = ZDX(4,2);
   pars(233) = ZDX(4,3);
   pars(234) = ZDX(4,4);
   pars(235) = ZDX(4,5);
   pars(236) = ZDX(5,0);
   pars(237) = ZDX(5,1);
   pars(238) = ZDX(5,2);
   pars(239) = ZDX(5,3);
   pars(240) = ZDX(5,4);
   pars(241) = ZDX(5,5);
   pars(242) = ZH(0,0);
   pars(243) = ZH(0,1);
   pars(244) = ZH(0,2);
   pars(245) = ZH(1,0);
   pars(246) = ZH(1,1);
   pars(247) = ZH(1,2);
   pars(248) = ZH(2,0);
   pars(249) = ZH(2,1);
   pars(250) = ZH(2,2);
   pars(251) = ZA(0,0);
   pars(252) = ZA(0,1);
   pars(253) = ZA(0,2);
   pars(254) = ZA(1,0);
   pars(255) = ZA(1,1);
   pars(256) = ZA(1,2);
   pars(257) = ZA(2,0);
   pars(258) = ZA(2,1);
   pars(259) = ZA(2,2);
   pars(260) = ZP(0,0);
   pars(261) = ZP(0,1);
   pars(262) = ZP(1,0);
   pars(263) = ZP(1,1);
   pars(264) = Re(ZN(0,0));
   pars(265) = Im(ZN(0,0));
   pars(266) = Re(ZN(0,1));
   pars(267) = Im(ZN(0,1));
   pars(268) = Re(ZN(0,2));
   pars(269) = Im(ZN(0,2));
   pars(270) = Re(ZN(0,3));
   pars(271) = Im(ZN(0,3));
   pars(272) = Re(ZN(0,4));
   pars(273) = Im(ZN(0,4));
   pars(274) = Re(ZN(0,5));
   pars(275) = Im(ZN(0,5));
   pars(276) = Re(ZN(1,0));
   pars(277) = Im(ZN(1,0));
   pars(278) = Re(ZN(1,1));
   pars(279) = Im(ZN(1,1));
   pars(280) = Re(ZN(1,2));
   pars(281) = Im(ZN(1,2));
   pars(282) = Re(ZN(1,3));
   pars(283) = Im(ZN(1,3));
   pars(284) = Re(ZN(1,4));
   pars(285) = Im(ZN(1,4));
   pars(286) = Re(ZN(1,5));
   pars(287) = Im(ZN(1,5));
   pars(288) = Re(ZN(2,0));
   pars(289) = Im(ZN(2,0));
   pars(290) = Re(ZN(2,1));
   pars(291) = Im(ZN(2,1));
   pars(292) = Re(ZN(2,2));
   pars(293) = Im(ZN(2,2));
   pars(294) = Re(ZN(2,3));
   pars(295) = Im(ZN(2,3));
   pars(296) = Re(ZN(2,4));
   pars(297) = Im(ZN(2,4));
   pars(298) = Re(ZN(2,5));
   pars(299) = Im(ZN(2,5));
   pars(300) = Re(ZN(3,0));
   pars(301) = Im(ZN(3,0));
   pars(302) = Re(ZN(3,1));
   pars(303) = Im(ZN(3,1));
   pars(304) = Re(ZN(3,2));
   pars(305) = Im(ZN(3,2));
   pars(306) = Re(ZN(3,3));
   pars(307) = Im(ZN(3,3));
   pars(308) = Re(ZN(3,4));
   pars(309) = Im(ZN(3,4));
   pars(310) = Re(ZN(3,5));
   pars(311) = Im(ZN(3,5));
   pars(312) = Re(ZN(4,0));
   pars(313) = Im(ZN(4,0));
   pars(314) = Re(ZN(4,1));
   pars(315) = Im(ZN(4,1));
   pars(316) = Re(ZN(4,2));
   pars(317) = Im(ZN(4,2));
   pars(318) = Re(ZN(4,3));
   pars(319) = Im(ZN(4,3));
   pars(320) = Re(ZN(4,4));
   pars(321) = Im(ZN(4,4));
   pars(322) = Re(ZN(4,5));
   pars(323) = Im(ZN(4,5));
   pars(324) = Re(ZN(5,0));
   pars(325) = Im(ZN(5,0));
   pars(326) = Re(ZN(5,1));
   pars(327) = Im(ZN(5,1));
   pars(328) = Re(ZN(5,2));
   pars(329) = Im(ZN(5,2));
   pars(330) = Re(ZN(5,3));
   pars(331) = Im(ZN(5,3));
   pars(332) = Re(ZN(5,4));
   pars(333) = Im(ZN(5,4));
   pars(334) = Re(ZN(5,5));
   pars(335) = Im(ZN(5,5));
   pars(336) = Re(UM(0,0));
   pars(337) = Im(UM(0,0));
   pars(338) = Re(UM(0,1));
   pars(339) = Im(UM(0,1));
   pars(340) = Re(UM(1,0));
   pars(341) = Im(UM(1,0));
   pars(342) = Re(UM(1,1));
   pars(343) = Im(UM(1,1));
   pars(344) = Re(UP(0,0));
   pars(345) = Im(UP(0,0));
   pars(346) = Re(UP(0,1));
   pars(347) = Im(UP(0,1));
   pars(348) = Re(UP(1,0));
   pars(349) = Im(UP(1,0));
   pars(350) = Re(UP(1,1));
   pars(351) = Im(UP(1,1));
   pars(352) = Re(ZEL(0,0));
   pars(353) = Im(ZEL(0,0));
   pars(354) = Re(ZEL(0,1));
   pars(355) = Im(ZEL(0,1));
   pars(356) = Re(ZEL(0,2));
   pars(357) = Im(ZEL(0,2));
   pars(358) = Re(ZEL(1,0));
   pars(359) = Im(ZEL(1,0));
   pars(360) = Re(ZEL(1,1));
   pars(361) = Im(ZEL(1,1));
   pars(362) = Re(ZEL(1,2));
   pars(363) = Im(ZEL(1,2));
   pars(364) = Re(ZEL(2,0));
   pars(365) = Im(ZEL(2,0));
   pars(366) = Re(ZEL(2,1));
   pars(367) = Im(ZEL(2,1));
   pars(368) = Re(ZEL(2,2));
   pars(369) = Im(ZEL(2,2));
   pars(370) = Re(ZER(0,0));
   pars(371) = Im(ZER(0,0));
   pars(372) = Re(ZER(0,1));
   pars(373) = Im(ZER(0,1));
   pars(374) = Re(ZER(0,2));
   pars(375) = Im(ZER(0,2));
   pars(376) = Re(ZER(1,0));
   pars(377) = Im(ZER(1,0));
   pars(378) = Re(ZER(1,1));
   pars(379) = Im(ZER(1,1));
   pars(380) = Re(ZER(1,2));
   pars(381) = Im(ZER(1,2));
   pars(382) = Re(ZER(2,0));
   pars(383) = Im(ZER(2,0));
   pars(384) = Re(ZER(2,1));
   pars(385) = Im(ZER(2,1));
   pars(386) = Re(ZER(2,2));
   pars(387) = Im(ZER(2,2));
   pars(388) = Re(ZDL(0,0));
   pars(389) = Im(ZDL(0,0));
   pars(390) = Re(ZDL(0,1));
   pars(391) = Im(ZDL(0,1));
   pars(392) = Re(ZDL(0,2));
   pars(393) = Im(ZDL(0,2));
   pars(394) = Re(ZDL(1,0));
   pars(395) = Im(ZDL(1,0));
   pars(396) = Re(ZDL(1,1));
   pars(397) = Im(ZDL(1,1));
   pars(398) = Re(ZDL(1,2));
   pars(399) = Im(ZDL(1,2));
   pars(400) = Re(ZDL(2,0));
   pars(401) = Im(ZDL(2,0));
   pars(402) = Re(ZDL(2,1));
   pars(403) = Im(ZDL(2,1));
   pars(404) = Re(ZDL(2,2));
   pars(405) = Im(ZDL(2,2));
   pars(406) = Re(ZDR(0,0));
   pars(407) = Im(ZDR(0,0));
   pars(408) = Re(ZDR(0,1));
   pars(409) = Im(ZDR(0,1));
   pars(410) = Re(ZDR(0,2));
   pars(411) = Im(ZDR(0,2));
   pars(412) = Re(ZDR(1,0));
   pars(413) = Im(ZDR(1,0));
   pars(414) = Re(ZDR(1,1));
   pars(415) = Im(ZDR(1,1));
   pars(416) = Re(ZDR(1,2));
   pars(417) = Im(ZDR(1,2));
   pars(418) = Re(ZDR(2,0));
   pars(419) = Im(ZDR(2,0));
   pars(420) = Re(ZDR(2,1));
   pars(421) = Im(ZDR(2,1));
   pars(422) = Re(ZDR(2,2));
   pars(423) = Im(ZDR(2,2));
   pars(424) = Re(ZUL(0,0));
   pars(425) = Im(ZUL(0,0));
   pars(426) = Re(ZUL(0,1));
   pars(427) = Im(ZUL(0,1));
   pars(428) = Re(ZUL(0,2));
   pars(429) = Im(ZUL(0,2));
   pars(430) = Re(ZUL(1,0));
   pars(431) = Im(ZUL(1,0));
   pars(432) = Re(ZUL(1,1));
   pars(433) = Im(ZUL(1,1));
   pars(434) = Re(ZUL(1,2));
   pars(435) = Im(ZUL(1,2));
   pars(436) = Re(ZUL(2,0));
   pars(437) = Im(ZUL(2,0));
   pars(438) = Re(ZUL(2,1));
   pars(439) = Im(ZUL(2,1));
   pars(440) = Re(ZUL(2,2));
   pars(441) = Im(ZUL(2,2));
   pars(442) = Re(ZUR(0,0));
   pars(443) = Im(ZUR(0,0));
   pars(444) = Re(ZUR(0,1));
   pars(445) = Im(ZUR(0,1));
   pars(446) = Re(ZUR(0,2));
   pars(447) = Im(ZUR(0,2));
   pars(448) = Re(ZUR(1,0));
   pars(449) = Im(ZUR(1,0));
   pars(450) = Re(ZUR(1,1));
   pars(451) = Im(ZUR(1,1));
   pars(452) = Re(ZUR(1,2));
   pars(453) = Im(ZUR(1,2));
   pars(454) = Re(ZUR(2,0));
   pars(455) = Im(ZUR(2,0));
   pars(456) = Re(ZUR(2,1));
   pars(457) = Im(ZUR(2,1));
   pars(458) = Re(ZUR(2,2));
   pars(459) = Im(ZUR(2,2));
   pars(460) = Re(ZDXL(0,0));
   pars(461) = Im(ZDXL(0,0));
   pars(462) = Re(ZDXL(0,1));
   pars(463) = Im(ZDXL(0,1));
   pars(464) = Re(ZDXL(0,2));
   pars(465) = Im(ZDXL(0,2));
   pars(466) = Re(ZDXL(1,0));
   pars(467) = Im(ZDXL(1,0));
   pars(468) = Re(ZDXL(1,1));
   pars(469) = Im(ZDXL(1,1));
   pars(470) = Re(ZDXL(1,2));
   pars(471) = Im(ZDXL(1,2));
   pars(472) = Re(ZDXL(2,0));
   pars(473) = Im(ZDXL(2,0));
   pars(474) = Re(ZDXL(2,1));
   pars(475) = Im(ZDXL(2,1));
   pars(476) = Re(ZDXL(2,2));
   pars(477) = Im(ZDXL(2,2));
   pars(478) = Re(ZDXR(0,0));
   pars(479) = Im(ZDXR(0,0));
   pars(480) = Re(ZDXR(0,1));
   pars(481) = Im(ZDXR(0,1));
   pars(482) = Re(ZDXR(0,2));
   pars(483) = Im(ZDXR(0,2));
   pars(484) = Re(ZDXR(1,0));
   pars(485) = Im(ZDXR(1,0));
   pars(486) = Re(ZDXR(1,1));
   pars(487) = Im(ZDXR(1,1));
   pars(488) = Re(ZDXR(1,2));
   pars(489) = Im(ZDXR(1,2));
   pars(490) = Re(ZDXR(2,0));
   pars(491) = Im(ZDXR(2,0));
   pars(492) = Re(ZDXR(2,1));
   pars(493) = Im(ZDXR(2,1));
   pars(494) = Re(ZDXR(2,2));
   pars(495) = Im(ZDXR(2,2));
   pars(496) = UHI0(0,0);
   pars(497) = UHI0(0,1);
   pars(498) = UHI0(0,2);
   pars(499) = UHI0(0,3);
   pars(500) = UHI0(1,0);
   pars(501) = UHI0(1,1);
   pars(502) = UHI0(1,2);
   pars(503) = UHI0(1,3);
   pars(504) = UHI0(2,0);
   pars(505) = UHI0(2,1);
   pars(506) = UHI0(2,2);
   pars(507) = UHI0(2,3);
   pars(508) = UHI0(3,0);
   pars(509) = UHI0(3,1);
   pars(510) = UHI0(3,2);
   pars(511) = UHI0(3,3);
   pars(512) = UHIp(0,0);
   pars(513) = UHIp(0,1);
   pars(514) = UHIp(0,2);
   pars(515) = UHIp(0,3);
   pars(516) = UHIp(1,0);
   pars(517) = UHIp(1,1);
   pars(518) = UHIp(1,2);
   pars(519) = UHIp(1,3);
   pars(520) = UHIp(2,0);
   pars(521) = UHIp(2,1);
   pars(522) = UHIp(2,2);
   pars(523) = UHIp(2,3);
   pars(524) = UHIp(3,0);
   pars(525) = UHIp(3,1);
   pars(526) = UHIp(3,2);
   pars(527) = UHIp(3,3);
   pars(528) = Re(ZMI(0,0));
   pars(529) = Im(ZMI(0,0));
   pars(530) = Re(ZMI(0,1));
   pars(531) = Im(ZMI(0,1));
   pars(532) = Re(ZMI(1,0));
   pars(533) = Im(ZMI(1,0));
   pars(534) = Re(ZMI(1,1));
   pars(535) = Im(ZMI(1,1));
   pars(536) = Re(ZPI(0,0));
   pars(537) = Im(ZPI(0,0));
   pars(538) = Re(ZPI(0,1));
   pars(539) = Im(ZPI(0,1));
   pars(540) = Re(ZPI(1,0));
   pars(541) = Im(ZPI(1,0));
   pars(542) = Re(ZPI(1,1));
   pars(543) = Im(ZPI(1,1));
   pars(544) = Re(ZNI(0,0));
   pars(545) = Im(ZNI(0,0));
   pars(546) = Re(ZNI(0,1));
   pars(547) = Im(ZNI(0,1));
   pars(548) = Re(ZNI(0,2));
   pars(549) = Im(ZNI(0,2));
   pars(550) = Re(ZNI(0,3));
   pars(551) = Im(ZNI(0,3));
   pars(552) = Re(ZNI(1,0));
   pars(553) = Im(ZNI(1,0));
   pars(554) = Re(ZNI(1,1));
   pars(555) = Im(ZNI(1,1));
   pars(556) = Re(ZNI(1,2));
   pars(557) = Im(ZNI(1,2));
   pars(558) = Re(ZNI(1,3));
   pars(559) = Im(ZNI(1,3));
   pars(560) = Re(ZNI(2,0));
   pars(561) = Im(ZNI(2,0));
   pars(562) = Re(ZNI(2,1));
   pars(563) = Im(ZNI(2,1));
   pars(564) = Re(ZNI(2,2));
   pars(565) = Im(ZNI(2,2));
   pars(566) = Re(ZNI(2,3));
   pars(567) = Im(ZNI(2,3));
   pars(568) = Re(ZNI(3,0));
   pars(569) = Im(ZNI(3,0));
   pars(570) = Re(ZNI(3,1));
   pars(571) = Im(ZNI(3,1));
   pars(572) = Re(ZNI(3,2));
   pars(573) = Im(ZNI(3,2));
   pars(574) = Re(ZNI(3,3));
   pars(575) = Im(ZNI(3,3));
   pars(576) = ZSSI(0,0);
   pars(577) = ZSSI(0,1);
   pars(578) = ZSSI(1,0);
   pars(579) = ZSSI(1,1);
   pars(580) = Re(ZFSI(0,0));
   pars(581) = Im(ZFSI(0,0));
   pars(582) = Re(ZFSI(0,1));
   pars(583) = Im(ZFSI(0,1));
   pars(584) = Re(ZFSI(1,0));
   pars(585) = Im(ZFSI(1,0));
   pars(586) = Re(ZFSI(1,1));
   pars(587) = Im(ZFSI(1,1));
   pars(588) = UHp0(0,0);
   pars(589) = UHp0(0,1);
   pars(590) = UHp0(1,0);
   pars(591) = UHp0(1,1);
   pars(592) = UHpp(0,0);
   pars(593) = UHpp(0,1);
   pars(594) = UHpp(1,0);
   pars(595) = UHpp(1,1);
   pars(596) = Re(ZNp(0,0));
   pars(597) = Im(ZNp(0,0));
   pars(598) = Re(ZNp(0,1));
   pars(599) = Im(ZNp(0,1));
   pars(600) = Re(ZNp(1,0));
   pars(601) = Im(ZNp(1,0));
   pars(602) = Re(ZNp(1,1));
   pars(603) = Im(ZNp(1,1));
   pars(604) = ZZ(0,0);
   pars(605) = ZZ(0,1);
   pars(606) = ZZ(0,2);
   pars(607) = ZZ(1,0);
   pars(608) = ZZ(1,1);
   pars(609) = ZZ(1,2);
   pars(610) = ZZ(2,0);
   pars(611) = ZZ(2,1);
   pars(612) = ZZ(2,2);


   return pars;
}

void E6SSMtower_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

   ZD(0,0) = pars(89);
   ZD(0,1) = pars(90);
   ZD(0,2) = pars(91);
   ZD(0,3) = pars(92);
   ZD(0,4) = pars(93);
   ZD(0,5) = pars(94);
   ZD(1,0) = pars(95);
   ZD(1,1) = pars(96);
   ZD(1,2) = pars(97);
   ZD(1,3) = pars(98);
   ZD(1,4) = pars(99);
   ZD(1,5) = pars(100);
   ZD(2,0) = pars(101);
   ZD(2,1) = pars(102);
   ZD(2,2) = pars(103);
   ZD(2,3) = pars(104);
   ZD(2,4) = pars(105);
   ZD(2,5) = pars(106);
   ZD(3,0) = pars(107);
   ZD(3,1) = pars(108);
   ZD(3,2) = pars(109);
   ZD(3,3) = pars(110);
   ZD(3,4) = pars(111);
   ZD(3,5) = pars(112);
   ZD(4,0) = pars(113);
   ZD(4,1) = pars(114);
   ZD(4,2) = pars(115);
   ZD(4,3) = pars(116);
   ZD(4,4) = pars(117);
   ZD(4,5) = pars(118);
   ZD(5,0) = pars(119);
   ZD(5,1) = pars(120);
   ZD(5,2) = pars(121);
   ZD(5,3) = pars(122);
   ZD(5,4) = pars(123);
   ZD(5,5) = pars(124);
   ZV(0,0) = pars(125);
   ZV(0,1) = pars(126);
   ZV(0,2) = pars(127);
   ZV(1,0) = pars(128);
   ZV(1,1) = pars(129);
   ZV(1,2) = pars(130);
   ZV(2,0) = pars(131);
   ZV(2,1) = pars(132);
   ZV(2,2) = pars(133);
   ZU(0,0) = pars(134);
   ZU(0,1) = pars(135);
   ZU(0,2) = pars(136);
   ZU(0,3) = pars(137);
   ZU(0,4) = pars(138);
   ZU(0,5) = pars(139);
   ZU(1,0) = pars(140);
   ZU(1,1) = pars(141);
   ZU(1,2) = pars(142);
   ZU(1,3) = pars(143);
   ZU(1,4) = pars(144);
   ZU(1,5) = pars(145);
   ZU(2,0) = pars(146);
   ZU(2,1) = pars(147);
   ZU(2,2) = pars(148);
   ZU(2,3) = pars(149);
   ZU(2,4) = pars(150);
   ZU(2,5) = pars(151);
   ZU(3,0) = pars(152);
   ZU(3,1) = pars(153);
   ZU(3,2) = pars(154);
   ZU(3,3) = pars(155);
   ZU(3,4) = pars(156);
   ZU(3,5) = pars(157);
   ZU(4,0) = pars(158);
   ZU(4,1) = pars(159);
   ZU(4,2) = pars(160);
   ZU(4,3) = pars(161);
   ZU(4,4) = pars(162);
   ZU(4,5) = pars(163);
   ZU(5,0) = pars(164);
   ZU(5,1) = pars(165);
   ZU(5,2) = pars(166);
   ZU(5,3) = pars(167);
   ZU(5,4) = pars(168);
   ZU(5,5) = pars(169);
   ZE(0,0) = pars(170);
   ZE(0,1) = pars(171);
   ZE(0,2) = pars(172);
   ZE(0,3) = pars(173);
   ZE(0,4) = pars(174);
   ZE(0,5) = pars(175);
   ZE(1,0) = pars(176);
   ZE(1,1) = pars(177);
   ZE(1,2) = pars(178);
   ZE(1,3) = pars(179);
   ZE(1,4) = pars(180);
   ZE(1,5) = pars(181);
   ZE(2,0) = pars(182);
   ZE(2,1) = pars(183);
   ZE(2,2) = pars(184);
   ZE(2,3) = pars(185);
   ZE(2,4) = pars(186);
   ZE(2,5) = pars(187);
   ZE(3,0) = pars(188);
   ZE(3,1) = pars(189);
   ZE(3,2) = pars(190);
   ZE(3,3) = pars(191);
   ZE(3,4) = pars(192);
   ZE(3,5) = pars(193);
   ZE(4,0) = pars(194);
   ZE(4,1) = pars(195);
   ZE(4,2) = pars(196);
   ZE(4,3) = pars(197);
   ZE(4,4) = pars(198);
   ZE(4,5) = pars(199);
   ZE(5,0) = pars(200);
   ZE(5,1) = pars(201);
   ZE(5,2) = pars(202);
   ZE(5,3) = pars(203);
   ZE(5,4) = pars(204);
   ZE(5,5) = pars(205);
   ZDX(0,0) = pars(206);
   ZDX(0,1) = pars(207);
   ZDX(0,2) = pars(208);
   ZDX(0,3) = pars(209);
   ZDX(0,4) = pars(210);
   ZDX(0,5) = pars(211);
   ZDX(1,0) = pars(212);
   ZDX(1,1) = pars(213);
   ZDX(1,2) = pars(214);
   ZDX(1,3) = pars(215);
   ZDX(1,4) = pars(216);
   ZDX(1,5) = pars(217);
   ZDX(2,0) = pars(218);
   ZDX(2,1) = pars(219);
   ZDX(2,2) = pars(220);
   ZDX(2,3) = pars(221);
   ZDX(2,4) = pars(222);
   ZDX(2,5) = pars(223);
   ZDX(3,0) = pars(224);
   ZDX(3,1) = pars(225);
   ZDX(3,2) = pars(226);
   ZDX(3,3) = pars(227);
   ZDX(3,4) = pars(228);
   ZDX(3,5) = pars(229);
   ZDX(4,0) = pars(230);
   ZDX(4,1) = pars(231);
   ZDX(4,2) = pars(232);
   ZDX(4,3) = pars(233);
   ZDX(4,4) = pars(234);
   ZDX(4,5) = pars(235);
   ZDX(5,0) = pars(236);
   ZDX(5,1) = pars(237);
   ZDX(5,2) = pars(238);
   ZDX(5,3) = pars(239);
   ZDX(5,4) = pars(240);
   ZDX(5,5) = pars(241);
   ZH(0,0) = pars(242);
   ZH(0,1) = pars(243);
   ZH(0,2) = pars(244);
   ZH(1,0) = pars(245);
   ZH(1,1) = pars(246);
   ZH(1,2) = pars(247);
   ZH(2,0) = pars(248);
   ZH(2,1) = pars(249);
   ZH(2,2) = pars(250);
   ZA(0,0) = pars(251);
   ZA(0,1) = pars(252);
   ZA(0,2) = pars(253);
   ZA(1,0) = pars(254);
   ZA(1,1) = pars(255);
   ZA(1,2) = pars(256);
   ZA(2,0) = pars(257);
   ZA(2,1) = pars(258);
   ZA(2,2) = pars(259);
   ZP(0,0) = pars(260);
   ZP(0,1) = pars(261);
   ZP(1,0) = pars(262);
   ZP(1,1) = pars(263);
   ZN(0,0) = std::complex<double>(pars(264), pars(265));
   ZN(0,1) = std::complex<double>(pars(266), pars(267));
   ZN(0,2) = std::complex<double>(pars(268), pars(269));
   ZN(0,3) = std::complex<double>(pars(270), pars(271));
   ZN(0,4) = std::complex<double>(pars(272), pars(273));
   ZN(0,5) = std::complex<double>(pars(274), pars(275));
   ZN(1,0) = std::complex<double>(pars(276), pars(277));
   ZN(1,1) = std::complex<double>(pars(278), pars(279));
   ZN(1,2) = std::complex<double>(pars(280), pars(281));
   ZN(1,3) = std::complex<double>(pars(282), pars(283));
   ZN(1,4) = std::complex<double>(pars(284), pars(285));
   ZN(1,5) = std::complex<double>(pars(286), pars(287));
   ZN(2,0) = std::complex<double>(pars(288), pars(289));
   ZN(2,1) = std::complex<double>(pars(290), pars(291));
   ZN(2,2) = std::complex<double>(pars(292), pars(293));
   ZN(2,3) = std::complex<double>(pars(294), pars(295));
   ZN(2,4) = std::complex<double>(pars(296), pars(297));
   ZN(2,5) = std::complex<double>(pars(298), pars(299));
   ZN(3,0) = std::complex<double>(pars(300), pars(301));
   ZN(3,1) = std::complex<double>(pars(302), pars(303));
   ZN(3,2) = std::complex<double>(pars(304), pars(305));
   ZN(3,3) = std::complex<double>(pars(306), pars(307));
   ZN(3,4) = std::complex<double>(pars(308), pars(309));
   ZN(3,5) = std::complex<double>(pars(310), pars(311));
   ZN(4,0) = std::complex<double>(pars(312), pars(313));
   ZN(4,1) = std::complex<double>(pars(314), pars(315));
   ZN(4,2) = std::complex<double>(pars(316), pars(317));
   ZN(4,3) = std::complex<double>(pars(318), pars(319));
   ZN(4,4) = std::complex<double>(pars(320), pars(321));
   ZN(4,5) = std::complex<double>(pars(322), pars(323));
   ZN(5,0) = std::complex<double>(pars(324), pars(325));
   ZN(5,1) = std::complex<double>(pars(326), pars(327));
   ZN(5,2) = std::complex<double>(pars(328), pars(329));
   ZN(5,3) = std::complex<double>(pars(330), pars(331));
   ZN(5,4) = std::complex<double>(pars(332), pars(333));
   ZN(5,5) = std::complex<double>(pars(334), pars(335));
   UM(0,0) = std::complex<double>(pars(336), pars(337));
   UM(0,1) = std::complex<double>(pars(338), pars(339));
   UM(1,0) = std::complex<double>(pars(340), pars(341));
   UM(1,1) = std::complex<double>(pars(342), pars(343));
   UP(0,0) = std::complex<double>(pars(344), pars(345));
   UP(0,1) = std::complex<double>(pars(346), pars(347));
   UP(1,0) = std::complex<double>(pars(348), pars(349));
   UP(1,1) = std::complex<double>(pars(350), pars(351));
   ZEL(0,0) = std::complex<double>(pars(352), pars(353));
   ZEL(0,1) = std::complex<double>(pars(354), pars(355));
   ZEL(0,2) = std::complex<double>(pars(356), pars(357));
   ZEL(1,0) = std::complex<double>(pars(358), pars(359));
   ZEL(1,1) = std::complex<double>(pars(360), pars(361));
   ZEL(1,2) = std::complex<double>(pars(362), pars(363));
   ZEL(2,0) = std::complex<double>(pars(364), pars(365));
   ZEL(2,1) = std::complex<double>(pars(366), pars(367));
   ZEL(2,2) = std::complex<double>(pars(368), pars(369));
   ZER(0,0) = std::complex<double>(pars(370), pars(371));
   ZER(0,1) = std::complex<double>(pars(372), pars(373));
   ZER(0,2) = std::complex<double>(pars(374), pars(375));
   ZER(1,0) = std::complex<double>(pars(376), pars(377));
   ZER(1,1) = std::complex<double>(pars(378), pars(379));
   ZER(1,2) = std::complex<double>(pars(380), pars(381));
   ZER(2,0) = std::complex<double>(pars(382), pars(383));
   ZER(2,1) = std::complex<double>(pars(384), pars(385));
   ZER(2,2) = std::complex<double>(pars(386), pars(387));
   ZDL(0,0) = std::complex<double>(pars(388), pars(389));
   ZDL(0,1) = std::complex<double>(pars(390), pars(391));
   ZDL(0,2) = std::complex<double>(pars(392), pars(393));
   ZDL(1,0) = std::complex<double>(pars(394), pars(395));
   ZDL(1,1) = std::complex<double>(pars(396), pars(397));
   ZDL(1,2) = std::complex<double>(pars(398), pars(399));
   ZDL(2,0) = std::complex<double>(pars(400), pars(401));
   ZDL(2,1) = std::complex<double>(pars(402), pars(403));
   ZDL(2,2) = std::complex<double>(pars(404), pars(405));
   ZDR(0,0) = std::complex<double>(pars(406), pars(407));
   ZDR(0,1) = std::complex<double>(pars(408), pars(409));
   ZDR(0,2) = std::complex<double>(pars(410), pars(411));
   ZDR(1,0) = std::complex<double>(pars(412), pars(413));
   ZDR(1,1) = std::complex<double>(pars(414), pars(415));
   ZDR(1,2) = std::complex<double>(pars(416), pars(417));
   ZDR(2,0) = std::complex<double>(pars(418), pars(419));
   ZDR(2,1) = std::complex<double>(pars(420), pars(421));
   ZDR(2,2) = std::complex<double>(pars(422), pars(423));
   ZUL(0,0) = std::complex<double>(pars(424), pars(425));
   ZUL(0,1) = std::complex<double>(pars(426), pars(427));
   ZUL(0,2) = std::complex<double>(pars(428), pars(429));
   ZUL(1,0) = std::complex<double>(pars(430), pars(431));
   ZUL(1,1) = std::complex<double>(pars(432), pars(433));
   ZUL(1,2) = std::complex<double>(pars(434), pars(435));
   ZUL(2,0) = std::complex<double>(pars(436), pars(437));
   ZUL(2,1) = std::complex<double>(pars(438), pars(439));
   ZUL(2,2) = std::complex<double>(pars(440), pars(441));
   ZUR(0,0) = std::complex<double>(pars(442), pars(443));
   ZUR(0,1) = std::complex<double>(pars(444), pars(445));
   ZUR(0,2) = std::complex<double>(pars(446), pars(447));
   ZUR(1,0) = std::complex<double>(pars(448), pars(449));
   ZUR(1,1) = std::complex<double>(pars(450), pars(451));
   ZUR(1,2) = std::complex<double>(pars(452), pars(453));
   ZUR(2,0) = std::complex<double>(pars(454), pars(455));
   ZUR(2,1) = std::complex<double>(pars(456), pars(457));
   ZUR(2,2) = std::complex<double>(pars(458), pars(459));
   ZDXL(0,0) = std::complex<double>(pars(460), pars(461));
   ZDXL(0,1) = std::complex<double>(pars(462), pars(463));
   ZDXL(0,2) = std::complex<double>(pars(464), pars(465));
   ZDXL(1,0) = std::complex<double>(pars(466), pars(467));
   ZDXL(1,1) = std::complex<double>(pars(468), pars(469));
   ZDXL(1,2) = std::complex<double>(pars(470), pars(471));
   ZDXL(2,0) = std::complex<double>(pars(472), pars(473));
   ZDXL(2,1) = std::complex<double>(pars(474), pars(475));
   ZDXL(2,2) = std::complex<double>(pars(476), pars(477));
   ZDXR(0,0) = std::complex<double>(pars(478), pars(479));
   ZDXR(0,1) = std::complex<double>(pars(480), pars(481));
   ZDXR(0,2) = std::complex<double>(pars(482), pars(483));
   ZDXR(1,0) = std::complex<double>(pars(484), pars(485));
   ZDXR(1,1) = std::complex<double>(pars(486), pars(487));
   ZDXR(1,2) = std::complex<double>(pars(488), pars(489));
   ZDXR(2,0) = std::complex<double>(pars(490), pars(491));
   ZDXR(2,1) = std::complex<double>(pars(492), pars(493));
   ZDXR(2,2) = std::complex<double>(pars(494), pars(495));
   UHI0(0,0) = pars(496);
   UHI0(0,1) = pars(497);
   UHI0(0,2) = pars(498);
   UHI0(0,3) = pars(499);
   UHI0(1,0) = pars(500);
   UHI0(1,1) = pars(501);
   UHI0(1,2) = pars(502);
   UHI0(1,3) = pars(503);
   UHI0(2,0) = pars(504);
   UHI0(2,1) = pars(505);
   UHI0(2,2) = pars(506);
   UHI0(2,3) = pars(507);
   UHI0(3,0) = pars(508);
   UHI0(3,1) = pars(509);
   UHI0(3,2) = pars(510);
   UHI0(3,3) = pars(511);
   UHIp(0,0) = pars(512);
   UHIp(0,1) = pars(513);
   UHIp(0,2) = pars(514);
   UHIp(0,3) = pars(515);
   UHIp(1,0) = pars(516);
   UHIp(1,1) = pars(517);
   UHIp(1,2) = pars(518);
   UHIp(1,3) = pars(519);
   UHIp(2,0) = pars(520);
   UHIp(2,1) = pars(521);
   UHIp(2,2) = pars(522);
   UHIp(2,3) = pars(523);
   UHIp(3,0) = pars(524);
   UHIp(3,1) = pars(525);
   UHIp(3,2) = pars(526);
   UHIp(3,3) = pars(527);
   ZMI(0,0) = std::complex<double>(pars(528), pars(529));
   ZMI(0,1) = std::complex<double>(pars(530), pars(531));
   ZMI(1,0) = std::complex<double>(pars(532), pars(533));
   ZMI(1,1) = std::complex<double>(pars(534), pars(535));
   ZPI(0,0) = std::complex<double>(pars(536), pars(537));
   ZPI(0,1) = std::complex<double>(pars(538), pars(539));
   ZPI(1,0) = std::complex<double>(pars(540), pars(541));
   ZPI(1,1) = std::complex<double>(pars(542), pars(543));
   ZNI(0,0) = std::complex<double>(pars(544), pars(545));
   ZNI(0,1) = std::complex<double>(pars(546), pars(547));
   ZNI(0,2) = std::complex<double>(pars(548), pars(549));
   ZNI(0,3) = std::complex<double>(pars(550), pars(551));
   ZNI(1,0) = std::complex<double>(pars(552), pars(553));
   ZNI(1,1) = std::complex<double>(pars(554), pars(555));
   ZNI(1,2) = std::complex<double>(pars(556), pars(557));
   ZNI(1,3) = std::complex<double>(pars(558), pars(559));
   ZNI(2,0) = std::complex<double>(pars(560), pars(561));
   ZNI(2,1) = std::complex<double>(pars(562), pars(563));
   ZNI(2,2) = std::complex<double>(pars(564), pars(565));
   ZNI(2,3) = std::complex<double>(pars(566), pars(567));
   ZNI(3,0) = std::complex<double>(pars(568), pars(569));
   ZNI(3,1) = std::complex<double>(pars(570), pars(571));
   ZNI(3,2) = std::complex<double>(pars(572), pars(573));
   ZNI(3,3) = std::complex<double>(pars(574), pars(575));
   ZSSI(0,0) = pars(576);
   ZSSI(0,1) = pars(577);
   ZSSI(1,0) = pars(578);
   ZSSI(1,1) = pars(579);
   ZFSI(0,0) = std::complex<double>(pars(580), pars(581));
   ZFSI(0,1) = std::complex<double>(pars(582), pars(583));
   ZFSI(1,0) = std::complex<double>(pars(584), pars(585));
   ZFSI(1,1) = std::complex<double>(pars(586), pars(587));
   UHp0(0,0) = pars(588);
   UHp0(0,1) = pars(589);
   UHp0(1,0) = pars(590);
   UHp0(1,1) = pars(591);
   UHpp(0,0) = pars(592);
   UHpp(0,1) = pars(593);
   UHpp(1,0) = pars(594);
   UHpp(1,1) = pars(595);
   ZNp(0,0) = std::complex<double>(pars(596), pars(597));
   ZNp(0,1) = std::complex<double>(pars(598), pars(599));
   ZNp(1,0) = std::complex<double>(pars(600), pars(601));
   ZNp(1,1) = std::complex<double>(pars(602), pars(603));
   ZZ(0,0) = pars(604);
   ZZ(0,1) = pars(605);
   ZZ(0,2) = pars(606);
   ZZ(1,0) = pars(607);
   ZZ(1,1) = pars(608);
   ZZ(1,2) = pars(609);
   ZZ(2,0) = pars(610);
   ZZ(2,1) = pars(611);
   ZZ(2,2) = pars(612);

}

Eigen::ArrayXd E6SSMtower_physical::get_masses() const
{
   Eigen::ArrayXd pars(89);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MChaP;
   pars(6) = MSd(0);
   pars(7) = MSd(1);
   pars(8) = MSd(2);
   pars(9) = MSd(3);
   pars(10) = MSd(4);
   pars(11) = MSd(5);
   pars(12) = MSv(0);
   pars(13) = MSv(1);
   pars(14) = MSv(2);
   pars(15) = MSu(0);
   pars(16) = MSu(1);
   pars(17) = MSu(2);
   pars(18) = MSu(3);
   pars(19) = MSu(4);
   pars(20) = MSu(5);
   pars(21) = MSe(0);
   pars(22) = MSe(1);
   pars(23) = MSe(2);
   pars(24) = MSe(3);
   pars(25) = MSe(4);
   pars(26) = MSe(5);
   pars(27) = MSDX(0);
   pars(28) = MSDX(1);
   pars(29) = MSDX(2);
   pars(30) = MSDX(3);
   pars(31) = MSDX(4);
   pars(32) = MSDX(5);
   pars(33) = Mhh(0);
   pars(34) = Mhh(1);
   pars(35) = Mhh(2);
   pars(36) = MAh(0);
   pars(37) = MAh(1);
   pars(38) = MAh(2);
   pars(39) = MHpm(0);
   pars(40) = MHpm(1);
   pars(41) = MChi(0);
   pars(42) = MChi(1);
   pars(43) = MChi(2);
   pars(44) = MChi(3);
   pars(45) = MChi(4);
   pars(46) = MChi(5);
   pars(47) = MCha(0);
   pars(48) = MCha(1);
   pars(49) = MFe(0);
   pars(50) = MFe(1);
   pars(51) = MFe(2);
   pars(52) = MFd(0);
   pars(53) = MFd(1);
   pars(54) = MFd(2);
   pars(55) = MFu(0);
   pars(56) = MFu(1);
   pars(57) = MFu(2);
   pars(58) = MFDX(0);
   pars(59) = MFDX(1);
   pars(60) = MFDX(2);
   pars(61) = MSHI0(0);
   pars(62) = MSHI0(1);
   pars(63) = MSHI0(2);
   pars(64) = MSHI0(3);
   pars(65) = MSHIp(0);
   pars(66) = MSHIp(1);
   pars(67) = MSHIp(2);
   pars(68) = MSHIp(3);
   pars(69) = MChaI(0);
   pars(70) = MChaI(1);
   pars(71) = MChiI(0);
   pars(72) = MChiI(1);
   pars(73) = MChiI(2);
   pars(74) = MChiI(3);
   pars(75) = MSSI0(0);
   pars(76) = MSSI0(1);
   pars(77) = MFSI(0);
   pars(78) = MFSI(1);
   pars(79) = MSHp0(0);
   pars(80) = MSHp0(1);
   pars(81) = MSHpp(0);
   pars(82) = MSHpp(1);
   pars(83) = MChiP(0);
   pars(84) = MChiP(1);
   pars(85) = MVWm;
   pars(86) = MVP;
   pars(87) = MVZ;
   pars(88) = MVZp;

   return pars;
}

void E6SSMtower_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MFv(0) = pars(2);
   MFv(1) = pars(3);
   MFv(2) = pars(4);
   MChaP = pars(5);
   MSd(0) = pars(6);
   MSd(1) = pars(7);
   MSd(2) = pars(8);
   MSd(3) = pars(9);
   MSd(4) = pars(10);
   MSd(5) = pars(11);
   MSv(0) = pars(12);
   MSv(1) = pars(13);
   MSv(2) = pars(14);
   MSu(0) = pars(15);
   MSu(1) = pars(16);
   MSu(2) = pars(17);
   MSu(3) = pars(18);
   MSu(4) = pars(19);
   MSu(5) = pars(20);
   MSe(0) = pars(21);
   MSe(1) = pars(22);
   MSe(2) = pars(23);
   MSe(3) = pars(24);
   MSe(4) = pars(25);
   MSe(5) = pars(26);
   MSDX(0) = pars(27);
   MSDX(1) = pars(28);
   MSDX(2) = pars(29);
   MSDX(3) = pars(30);
   MSDX(4) = pars(31);
   MSDX(5) = pars(32);
   Mhh(0) = pars(33);
   Mhh(1) = pars(34);
   Mhh(2) = pars(35);
   MAh(0) = pars(36);
   MAh(1) = pars(37);
   MAh(2) = pars(38);
   MHpm(0) = pars(39);
   MHpm(1) = pars(40);
   MChi(0) = pars(41);
   MChi(1) = pars(42);
   MChi(2) = pars(43);
   MChi(3) = pars(44);
   MChi(4) = pars(45);
   MChi(5) = pars(46);
   MCha(0) = pars(47);
   MCha(1) = pars(48);
   MFe(0) = pars(49);
   MFe(1) = pars(50);
   MFe(2) = pars(51);
   MFd(0) = pars(52);
   MFd(1) = pars(53);
   MFd(2) = pars(54);
   MFu(0) = pars(55);
   MFu(1) = pars(56);
   MFu(2) = pars(57);
   MFDX(0) = pars(58);
   MFDX(1) = pars(59);
   MFDX(2) = pars(60);
   MSHI0(0) = pars(61);
   MSHI0(1) = pars(62);
   MSHI0(2) = pars(63);
   MSHI0(3) = pars(64);
   MSHIp(0) = pars(65);
   MSHIp(1) = pars(66);
   MSHIp(2) = pars(67);
   MSHIp(3) = pars(68);
   MChaI(0) = pars(69);
   MChaI(1) = pars(70);
   MChiI(0) = pars(71);
   MChiI(1) = pars(72);
   MChiI(2) = pars(73);
   MChiI(3) = pars(74);
   MSSI0(0) = pars(75);
   MSSI0(1) = pars(76);
   MFSI(0) = pars(77);
   MFSI(1) = pars(78);
   MSHp0(0) = pars(79);
   MSHp0(1) = pars(80);
   MSHpp(0) = pars(81);
   MSHpp(1) = pars(82);
   MChiP(0) = pars(83);
   MChiP(1) = pars(84);
   MVWm = pars(85);
   MVP = pars(86);
   MVZ = pars(87);
   MVZp = pars(88);

}

void E6SSMtower_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MChaP = " << MChaP << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "MSDX = " << MSDX.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFDX = " << MFDX.transpose() << '\n';
   ostr << "MSHI0 = " << MSHI0.transpose() << '\n';
   ostr << "MSHIp = " << MSHIp.transpose() << '\n';
   ostr << "MChaI = " << MChaI.transpose() << '\n';
   ostr << "MChiI = " << MChiI.transpose() << '\n';
   ostr << "MSSI0 = " << MSSI0.transpose() << '\n';
   ostr << "MFSI = " << MFSI.transpose() << '\n';
   ostr << "MSHp0 = " << MSHp0.transpose() << '\n';
   ostr << "MSHpp = " << MSHpp.transpose() << '\n';
   ostr << "MChiP = " << MChiP.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';

   ostr << "----------------------------------------\n"
           "pole mass mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZV = " << ZV << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZDX = " << ZDX << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UM = " << UM << '\n';
   ostr << "UP = " << UP << '\n';
   ostr << "ZEL = " << ZEL << '\n';
   ostr << "ZER = " << ZER << '\n';
   ostr << "ZDL = " << ZDL << '\n';
   ostr << "ZDR = " << ZDR << '\n';
   ostr << "ZUL = " << ZUL << '\n';
   ostr << "ZUR = " << ZUR << '\n';
   ostr << "ZDXL = " << ZDXL << '\n';
   ostr << "ZDXR = " << ZDXR << '\n';
   ostr << "UHI0 = " << UHI0 << '\n';
   ostr << "UHIp = " << UHIp << '\n';
   ostr << "ZMI = " << ZMI << '\n';
   ostr << "ZPI = " << ZPI << '\n';
   ostr << "ZNI = " << ZNI << '\n';
   ostr << "ZSSI = " << ZSSI << '\n';
   ostr << "ZFSI = " << ZFSI << '\n';
   ostr << "UHp0 = " << UHp0 << '\n';
   ostr << "UHpp = " << UHpp << '\n';
   ostr << "ZNp = " << ZNp << '\n';
   ostr << "ZZ = " << ZZ << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const E6SSMtower_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
