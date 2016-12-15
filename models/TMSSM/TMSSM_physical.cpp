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

// File generated at Thu 15 Dec 2016 12:47:25

#include "TMSSM_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

TMSSM_physical::TMSSM_physical()
   :
    MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MSd(Eigen::Array<
       double,6,1>::Zero()), MSv(Eigen::Array<double,3,1>::Zero()), MSu(
       Eigen::Array<double,6,1>::Zero()), MSe(Eigen::Array<double,6,1>::Zero()),
       Mhh(Eigen::Array<double,3,1>::Zero()), MAh(Eigen::Array<double,3,1>::Zero()
       ), MHpm(Eigen::Array<double,4,1>::Zero()), MChi(Eigen::Array<double,5,1>
       ::Zero()), MCha(Eigen::Array<double,3,1>::Zero()), MFe(Eigen::Array<double,
       3,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<
       double,3,1>::Zero()), MVWm(0), MVP(0), MVZ(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,3,3>::Zero()), ZA(Eigen::Matrix<double,3,
      3>::Zero()), ZP(Eigen::Matrix<double,4,4>::Zero()), ZN(Eigen::Matrix<
      std::complex<double>,5,5>::Zero()), UM(Eigen::Matrix<std::complex<double>,3,
      3>::Zero()), UP(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZEL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZER(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZDL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZDR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())

{
}

void TMSSM_physical::clear()
{
   MVG = 0.;
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,4,1>::Zero();
   ZP = Eigen::Matrix<double,4,4>::Zero();
   MChi = Eigen::Matrix<double,5,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,5,5>::Zero();
   MCha = Eigen::Matrix<double,3,1>::Zero();
   UM = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   UP = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   ZEL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZER = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFd = Eigen::Matrix<double,3,1>::Zero();
   ZDL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZDR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   ZUL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZUR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MVWm = 0.;
   MVP = 0.;
   MVZ = 0.;

}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void TMSSM_physical::convert_to_hk()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void TMSSM_physical::convert_to_slha()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

Eigen::ArrayXd TMSSM_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(405);

   pars(56) = ZD(0,0);
   pars(57) = ZD(0,1);
   pars(58) = ZD(0,2);
   pars(59) = ZD(0,3);
   pars(60) = ZD(0,4);
   pars(61) = ZD(0,5);
   pars(62) = ZD(1,0);
   pars(63) = ZD(1,1);
   pars(64) = ZD(1,2);
   pars(65) = ZD(1,3);
   pars(66) = ZD(1,4);
   pars(67) = ZD(1,5);
   pars(68) = ZD(2,0);
   pars(69) = ZD(2,1);
   pars(70) = ZD(2,2);
   pars(71) = ZD(2,3);
   pars(72) = ZD(2,4);
   pars(73) = ZD(2,5);
   pars(74) = ZD(3,0);
   pars(75) = ZD(3,1);
   pars(76) = ZD(3,2);
   pars(77) = ZD(3,3);
   pars(78) = ZD(3,4);
   pars(79) = ZD(3,5);
   pars(80) = ZD(4,0);
   pars(81) = ZD(4,1);
   pars(82) = ZD(4,2);
   pars(83) = ZD(4,3);
   pars(84) = ZD(4,4);
   pars(85) = ZD(4,5);
   pars(86) = ZD(5,0);
   pars(87) = ZD(5,1);
   pars(88) = ZD(5,2);
   pars(89) = ZD(5,3);
   pars(90) = ZD(5,4);
   pars(91) = ZD(5,5);
   pars(92) = ZV(0,0);
   pars(93) = ZV(0,1);
   pars(94) = ZV(0,2);
   pars(95) = ZV(1,0);
   pars(96) = ZV(1,1);
   pars(97) = ZV(1,2);
   pars(98) = ZV(2,0);
   pars(99) = ZV(2,1);
   pars(100) = ZV(2,2);
   pars(101) = ZU(0,0);
   pars(102) = ZU(0,1);
   pars(103) = ZU(0,2);
   pars(104) = ZU(0,3);
   pars(105) = ZU(0,4);
   pars(106) = ZU(0,5);
   pars(107) = ZU(1,0);
   pars(108) = ZU(1,1);
   pars(109) = ZU(1,2);
   pars(110) = ZU(1,3);
   pars(111) = ZU(1,4);
   pars(112) = ZU(1,5);
   pars(113) = ZU(2,0);
   pars(114) = ZU(2,1);
   pars(115) = ZU(2,2);
   pars(116) = ZU(2,3);
   pars(117) = ZU(2,4);
   pars(118) = ZU(2,5);
   pars(119) = ZU(3,0);
   pars(120) = ZU(3,1);
   pars(121) = ZU(3,2);
   pars(122) = ZU(3,3);
   pars(123) = ZU(3,4);
   pars(124) = ZU(3,5);
   pars(125) = ZU(4,0);
   pars(126) = ZU(4,1);
   pars(127) = ZU(4,2);
   pars(128) = ZU(4,3);
   pars(129) = ZU(4,4);
   pars(130) = ZU(4,5);
   pars(131) = ZU(5,0);
   pars(132) = ZU(5,1);
   pars(133) = ZU(5,2);
   pars(134) = ZU(5,3);
   pars(135) = ZU(5,4);
   pars(136) = ZU(5,5);
   pars(137) = ZE(0,0);
   pars(138) = ZE(0,1);
   pars(139) = ZE(0,2);
   pars(140) = ZE(0,3);
   pars(141) = ZE(0,4);
   pars(142) = ZE(0,5);
   pars(143) = ZE(1,0);
   pars(144) = ZE(1,1);
   pars(145) = ZE(1,2);
   pars(146) = ZE(1,3);
   pars(147) = ZE(1,4);
   pars(148) = ZE(1,5);
   pars(149) = ZE(2,0);
   pars(150) = ZE(2,1);
   pars(151) = ZE(2,2);
   pars(152) = ZE(2,3);
   pars(153) = ZE(2,4);
   pars(154) = ZE(2,5);
   pars(155) = ZE(3,0);
   pars(156) = ZE(3,1);
   pars(157) = ZE(3,2);
   pars(158) = ZE(3,3);
   pars(159) = ZE(3,4);
   pars(160) = ZE(3,5);
   pars(161) = ZE(4,0);
   pars(162) = ZE(4,1);
   pars(163) = ZE(4,2);
   pars(164) = ZE(4,3);
   pars(165) = ZE(4,4);
   pars(166) = ZE(4,5);
   pars(167) = ZE(5,0);
   pars(168) = ZE(5,1);
   pars(169) = ZE(5,2);
   pars(170) = ZE(5,3);
   pars(171) = ZE(5,4);
   pars(172) = ZE(5,5);
   pars(173) = ZH(0,0);
   pars(174) = ZH(0,1);
   pars(175) = ZH(0,2);
   pars(176) = ZH(1,0);
   pars(177) = ZH(1,1);
   pars(178) = ZH(1,2);
   pars(179) = ZH(2,0);
   pars(180) = ZH(2,1);
   pars(181) = ZH(2,2);
   pars(182) = ZA(0,0);
   pars(183) = ZA(0,1);
   pars(184) = ZA(0,2);
   pars(185) = ZA(1,0);
   pars(186) = ZA(1,1);
   pars(187) = ZA(1,2);
   pars(188) = ZA(2,0);
   pars(189) = ZA(2,1);
   pars(190) = ZA(2,2);
   pars(191) = ZP(0,0);
   pars(192) = ZP(0,1);
   pars(193) = ZP(0,2);
   pars(194) = ZP(0,3);
   pars(195) = ZP(1,0);
   pars(196) = ZP(1,1);
   pars(197) = ZP(1,2);
   pars(198) = ZP(1,3);
   pars(199) = ZP(2,0);
   pars(200) = ZP(2,1);
   pars(201) = ZP(2,2);
   pars(202) = ZP(2,3);
   pars(203) = ZP(3,0);
   pars(204) = ZP(3,1);
   pars(205) = ZP(3,2);
   pars(206) = ZP(3,3);
   pars(207) = Re(ZN(0,0));
   pars(208) = Im(ZN(0,0));
   pars(209) = Re(ZN(0,1));
   pars(210) = Im(ZN(0,1));
   pars(211) = Re(ZN(0,2));
   pars(212) = Im(ZN(0,2));
   pars(213) = Re(ZN(0,3));
   pars(214) = Im(ZN(0,3));
   pars(215) = Re(ZN(0,4));
   pars(216) = Im(ZN(0,4));
   pars(217) = Re(ZN(1,0));
   pars(218) = Im(ZN(1,0));
   pars(219) = Re(ZN(1,1));
   pars(220) = Im(ZN(1,1));
   pars(221) = Re(ZN(1,2));
   pars(222) = Im(ZN(1,2));
   pars(223) = Re(ZN(1,3));
   pars(224) = Im(ZN(1,3));
   pars(225) = Re(ZN(1,4));
   pars(226) = Im(ZN(1,4));
   pars(227) = Re(ZN(2,0));
   pars(228) = Im(ZN(2,0));
   pars(229) = Re(ZN(2,1));
   pars(230) = Im(ZN(2,1));
   pars(231) = Re(ZN(2,2));
   pars(232) = Im(ZN(2,2));
   pars(233) = Re(ZN(2,3));
   pars(234) = Im(ZN(2,3));
   pars(235) = Re(ZN(2,4));
   pars(236) = Im(ZN(2,4));
   pars(237) = Re(ZN(3,0));
   pars(238) = Im(ZN(3,0));
   pars(239) = Re(ZN(3,1));
   pars(240) = Im(ZN(3,1));
   pars(241) = Re(ZN(3,2));
   pars(242) = Im(ZN(3,2));
   pars(243) = Re(ZN(3,3));
   pars(244) = Im(ZN(3,3));
   pars(245) = Re(ZN(3,4));
   pars(246) = Im(ZN(3,4));
   pars(247) = Re(ZN(4,0));
   pars(248) = Im(ZN(4,0));
   pars(249) = Re(ZN(4,1));
   pars(250) = Im(ZN(4,1));
   pars(251) = Re(ZN(4,2));
   pars(252) = Im(ZN(4,2));
   pars(253) = Re(ZN(4,3));
   pars(254) = Im(ZN(4,3));
   pars(255) = Re(ZN(4,4));
   pars(256) = Im(ZN(4,4));
   pars(257) = Re(UM(0,0));
   pars(258) = Im(UM(0,0));
   pars(259) = Re(UM(0,1));
   pars(260) = Im(UM(0,1));
   pars(261) = Re(UM(0,2));
   pars(262) = Im(UM(0,2));
   pars(263) = Re(UM(1,0));
   pars(264) = Im(UM(1,0));
   pars(265) = Re(UM(1,1));
   pars(266) = Im(UM(1,1));
   pars(267) = Re(UM(1,2));
   pars(268) = Im(UM(1,2));
   pars(269) = Re(UM(2,0));
   pars(270) = Im(UM(2,0));
   pars(271) = Re(UM(2,1));
   pars(272) = Im(UM(2,1));
   pars(273) = Re(UM(2,2));
   pars(274) = Im(UM(2,2));
   pars(275) = Re(UP(0,0));
   pars(276) = Im(UP(0,0));
   pars(277) = Re(UP(0,1));
   pars(278) = Im(UP(0,1));
   pars(279) = Re(UP(0,2));
   pars(280) = Im(UP(0,2));
   pars(281) = Re(UP(1,0));
   pars(282) = Im(UP(1,0));
   pars(283) = Re(UP(1,1));
   pars(284) = Im(UP(1,1));
   pars(285) = Re(UP(1,2));
   pars(286) = Im(UP(1,2));
   pars(287) = Re(UP(2,0));
   pars(288) = Im(UP(2,0));
   pars(289) = Re(UP(2,1));
   pars(290) = Im(UP(2,1));
   pars(291) = Re(UP(2,2));
   pars(292) = Im(UP(2,2));
   pars(293) = Re(ZEL(0,0));
   pars(294) = Im(ZEL(0,0));
   pars(295) = Re(ZEL(0,1));
   pars(296) = Im(ZEL(0,1));
   pars(297) = Re(ZEL(0,2));
   pars(298) = Im(ZEL(0,2));
   pars(299) = Re(ZEL(1,0));
   pars(300) = Im(ZEL(1,0));
   pars(301) = Re(ZEL(1,1));
   pars(302) = Im(ZEL(1,1));
   pars(303) = Re(ZEL(1,2));
   pars(304) = Im(ZEL(1,2));
   pars(305) = Re(ZEL(2,0));
   pars(306) = Im(ZEL(2,0));
   pars(307) = Re(ZEL(2,1));
   pars(308) = Im(ZEL(2,1));
   pars(309) = Re(ZEL(2,2));
   pars(310) = Im(ZEL(2,2));
   pars(311) = Re(ZER(0,0));
   pars(312) = Im(ZER(0,0));
   pars(313) = Re(ZER(0,1));
   pars(314) = Im(ZER(0,1));
   pars(315) = Re(ZER(0,2));
   pars(316) = Im(ZER(0,2));
   pars(317) = Re(ZER(1,0));
   pars(318) = Im(ZER(1,0));
   pars(319) = Re(ZER(1,1));
   pars(320) = Im(ZER(1,1));
   pars(321) = Re(ZER(1,2));
   pars(322) = Im(ZER(1,2));
   pars(323) = Re(ZER(2,0));
   pars(324) = Im(ZER(2,0));
   pars(325) = Re(ZER(2,1));
   pars(326) = Im(ZER(2,1));
   pars(327) = Re(ZER(2,2));
   pars(328) = Im(ZER(2,2));
   pars(329) = Re(ZDL(0,0));
   pars(330) = Im(ZDL(0,0));
   pars(331) = Re(ZDL(0,1));
   pars(332) = Im(ZDL(0,1));
   pars(333) = Re(ZDL(0,2));
   pars(334) = Im(ZDL(0,2));
   pars(335) = Re(ZDL(1,0));
   pars(336) = Im(ZDL(1,0));
   pars(337) = Re(ZDL(1,1));
   pars(338) = Im(ZDL(1,1));
   pars(339) = Re(ZDL(1,2));
   pars(340) = Im(ZDL(1,2));
   pars(341) = Re(ZDL(2,0));
   pars(342) = Im(ZDL(2,0));
   pars(343) = Re(ZDL(2,1));
   pars(344) = Im(ZDL(2,1));
   pars(345) = Re(ZDL(2,2));
   pars(346) = Im(ZDL(2,2));
   pars(347) = Re(ZDR(0,0));
   pars(348) = Im(ZDR(0,0));
   pars(349) = Re(ZDR(0,1));
   pars(350) = Im(ZDR(0,1));
   pars(351) = Re(ZDR(0,2));
   pars(352) = Im(ZDR(0,2));
   pars(353) = Re(ZDR(1,0));
   pars(354) = Im(ZDR(1,0));
   pars(355) = Re(ZDR(1,1));
   pars(356) = Im(ZDR(1,1));
   pars(357) = Re(ZDR(1,2));
   pars(358) = Im(ZDR(1,2));
   pars(359) = Re(ZDR(2,0));
   pars(360) = Im(ZDR(2,0));
   pars(361) = Re(ZDR(2,1));
   pars(362) = Im(ZDR(2,1));
   pars(363) = Re(ZDR(2,2));
   pars(364) = Im(ZDR(2,2));
   pars(365) = Re(ZUL(0,0));
   pars(366) = Im(ZUL(0,0));
   pars(367) = Re(ZUL(0,1));
   pars(368) = Im(ZUL(0,1));
   pars(369) = Re(ZUL(0,2));
   pars(370) = Im(ZUL(0,2));
   pars(371) = Re(ZUL(1,0));
   pars(372) = Im(ZUL(1,0));
   pars(373) = Re(ZUL(1,1));
   pars(374) = Im(ZUL(1,1));
   pars(375) = Re(ZUL(1,2));
   pars(376) = Im(ZUL(1,2));
   pars(377) = Re(ZUL(2,0));
   pars(378) = Im(ZUL(2,0));
   pars(379) = Re(ZUL(2,1));
   pars(380) = Im(ZUL(2,1));
   pars(381) = Re(ZUL(2,2));
   pars(382) = Im(ZUL(2,2));
   pars(383) = Re(ZUR(0,0));
   pars(384) = Im(ZUR(0,0));
   pars(385) = Re(ZUR(0,1));
   pars(386) = Im(ZUR(0,1));
   pars(387) = Re(ZUR(0,2));
   pars(388) = Im(ZUR(0,2));
   pars(389) = Re(ZUR(1,0));
   pars(390) = Im(ZUR(1,0));
   pars(391) = Re(ZUR(1,1));
   pars(392) = Im(ZUR(1,1));
   pars(393) = Re(ZUR(1,2));
   pars(394) = Im(ZUR(1,2));
   pars(395) = Re(ZUR(2,0));
   pars(396) = Im(ZUR(2,0));
   pars(397) = Re(ZUR(2,1));
   pars(398) = Im(ZUR(2,1));
   pars(399) = Re(ZUR(2,2));
   pars(400) = Im(ZUR(2,2));
   pars(401) = ZZ(0,0);
   pars(402) = ZZ(0,1);
   pars(403) = ZZ(1,0);
   pars(404) = ZZ(1,1);


   return pars;
}

void TMSSM_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

   ZD(0,0) = pars(56);
   ZD(0,1) = pars(57);
   ZD(0,2) = pars(58);
   ZD(0,3) = pars(59);
   ZD(0,4) = pars(60);
   ZD(0,5) = pars(61);
   ZD(1,0) = pars(62);
   ZD(1,1) = pars(63);
   ZD(1,2) = pars(64);
   ZD(1,3) = pars(65);
   ZD(1,4) = pars(66);
   ZD(1,5) = pars(67);
   ZD(2,0) = pars(68);
   ZD(2,1) = pars(69);
   ZD(2,2) = pars(70);
   ZD(2,3) = pars(71);
   ZD(2,4) = pars(72);
   ZD(2,5) = pars(73);
   ZD(3,0) = pars(74);
   ZD(3,1) = pars(75);
   ZD(3,2) = pars(76);
   ZD(3,3) = pars(77);
   ZD(3,4) = pars(78);
   ZD(3,5) = pars(79);
   ZD(4,0) = pars(80);
   ZD(4,1) = pars(81);
   ZD(4,2) = pars(82);
   ZD(4,3) = pars(83);
   ZD(4,4) = pars(84);
   ZD(4,5) = pars(85);
   ZD(5,0) = pars(86);
   ZD(5,1) = pars(87);
   ZD(5,2) = pars(88);
   ZD(5,3) = pars(89);
   ZD(5,4) = pars(90);
   ZD(5,5) = pars(91);
   ZV(0,0) = pars(92);
   ZV(0,1) = pars(93);
   ZV(0,2) = pars(94);
   ZV(1,0) = pars(95);
   ZV(1,1) = pars(96);
   ZV(1,2) = pars(97);
   ZV(2,0) = pars(98);
   ZV(2,1) = pars(99);
   ZV(2,2) = pars(100);
   ZU(0,0) = pars(101);
   ZU(0,1) = pars(102);
   ZU(0,2) = pars(103);
   ZU(0,3) = pars(104);
   ZU(0,4) = pars(105);
   ZU(0,5) = pars(106);
   ZU(1,0) = pars(107);
   ZU(1,1) = pars(108);
   ZU(1,2) = pars(109);
   ZU(1,3) = pars(110);
   ZU(1,4) = pars(111);
   ZU(1,5) = pars(112);
   ZU(2,0) = pars(113);
   ZU(2,1) = pars(114);
   ZU(2,2) = pars(115);
   ZU(2,3) = pars(116);
   ZU(2,4) = pars(117);
   ZU(2,5) = pars(118);
   ZU(3,0) = pars(119);
   ZU(3,1) = pars(120);
   ZU(3,2) = pars(121);
   ZU(3,3) = pars(122);
   ZU(3,4) = pars(123);
   ZU(3,5) = pars(124);
   ZU(4,0) = pars(125);
   ZU(4,1) = pars(126);
   ZU(4,2) = pars(127);
   ZU(4,3) = pars(128);
   ZU(4,4) = pars(129);
   ZU(4,5) = pars(130);
   ZU(5,0) = pars(131);
   ZU(5,1) = pars(132);
   ZU(5,2) = pars(133);
   ZU(5,3) = pars(134);
   ZU(5,4) = pars(135);
   ZU(5,5) = pars(136);
   ZE(0,0) = pars(137);
   ZE(0,1) = pars(138);
   ZE(0,2) = pars(139);
   ZE(0,3) = pars(140);
   ZE(0,4) = pars(141);
   ZE(0,5) = pars(142);
   ZE(1,0) = pars(143);
   ZE(1,1) = pars(144);
   ZE(1,2) = pars(145);
   ZE(1,3) = pars(146);
   ZE(1,4) = pars(147);
   ZE(1,5) = pars(148);
   ZE(2,0) = pars(149);
   ZE(2,1) = pars(150);
   ZE(2,2) = pars(151);
   ZE(2,3) = pars(152);
   ZE(2,4) = pars(153);
   ZE(2,5) = pars(154);
   ZE(3,0) = pars(155);
   ZE(3,1) = pars(156);
   ZE(3,2) = pars(157);
   ZE(3,3) = pars(158);
   ZE(3,4) = pars(159);
   ZE(3,5) = pars(160);
   ZE(4,0) = pars(161);
   ZE(4,1) = pars(162);
   ZE(4,2) = pars(163);
   ZE(4,3) = pars(164);
   ZE(4,4) = pars(165);
   ZE(4,5) = pars(166);
   ZE(5,0) = pars(167);
   ZE(5,1) = pars(168);
   ZE(5,2) = pars(169);
   ZE(5,3) = pars(170);
   ZE(5,4) = pars(171);
   ZE(5,5) = pars(172);
   ZH(0,0) = pars(173);
   ZH(0,1) = pars(174);
   ZH(0,2) = pars(175);
   ZH(1,0) = pars(176);
   ZH(1,1) = pars(177);
   ZH(1,2) = pars(178);
   ZH(2,0) = pars(179);
   ZH(2,1) = pars(180);
   ZH(2,2) = pars(181);
   ZA(0,0) = pars(182);
   ZA(0,1) = pars(183);
   ZA(0,2) = pars(184);
   ZA(1,0) = pars(185);
   ZA(1,1) = pars(186);
   ZA(1,2) = pars(187);
   ZA(2,0) = pars(188);
   ZA(2,1) = pars(189);
   ZA(2,2) = pars(190);
   ZP(0,0) = pars(191);
   ZP(0,1) = pars(192);
   ZP(0,2) = pars(193);
   ZP(0,3) = pars(194);
   ZP(1,0) = pars(195);
   ZP(1,1) = pars(196);
   ZP(1,2) = pars(197);
   ZP(1,3) = pars(198);
   ZP(2,0) = pars(199);
   ZP(2,1) = pars(200);
   ZP(2,2) = pars(201);
   ZP(2,3) = pars(202);
   ZP(3,0) = pars(203);
   ZP(3,1) = pars(204);
   ZP(3,2) = pars(205);
   ZP(3,3) = pars(206);
   ZN(0,0) = std::complex<double>(pars(207), pars(208));
   ZN(0,1) = std::complex<double>(pars(209), pars(210));
   ZN(0,2) = std::complex<double>(pars(211), pars(212));
   ZN(0,3) = std::complex<double>(pars(213), pars(214));
   ZN(0,4) = std::complex<double>(pars(215), pars(216));
   ZN(1,0) = std::complex<double>(pars(217), pars(218));
   ZN(1,1) = std::complex<double>(pars(219), pars(220));
   ZN(1,2) = std::complex<double>(pars(221), pars(222));
   ZN(1,3) = std::complex<double>(pars(223), pars(224));
   ZN(1,4) = std::complex<double>(pars(225), pars(226));
   ZN(2,0) = std::complex<double>(pars(227), pars(228));
   ZN(2,1) = std::complex<double>(pars(229), pars(230));
   ZN(2,2) = std::complex<double>(pars(231), pars(232));
   ZN(2,3) = std::complex<double>(pars(233), pars(234));
   ZN(2,4) = std::complex<double>(pars(235), pars(236));
   ZN(3,0) = std::complex<double>(pars(237), pars(238));
   ZN(3,1) = std::complex<double>(pars(239), pars(240));
   ZN(3,2) = std::complex<double>(pars(241), pars(242));
   ZN(3,3) = std::complex<double>(pars(243), pars(244));
   ZN(3,4) = std::complex<double>(pars(245), pars(246));
   ZN(4,0) = std::complex<double>(pars(247), pars(248));
   ZN(4,1) = std::complex<double>(pars(249), pars(250));
   ZN(4,2) = std::complex<double>(pars(251), pars(252));
   ZN(4,3) = std::complex<double>(pars(253), pars(254));
   ZN(4,4) = std::complex<double>(pars(255), pars(256));
   UM(0,0) = std::complex<double>(pars(257), pars(258));
   UM(0,1) = std::complex<double>(pars(259), pars(260));
   UM(0,2) = std::complex<double>(pars(261), pars(262));
   UM(1,0) = std::complex<double>(pars(263), pars(264));
   UM(1,1) = std::complex<double>(pars(265), pars(266));
   UM(1,2) = std::complex<double>(pars(267), pars(268));
   UM(2,0) = std::complex<double>(pars(269), pars(270));
   UM(2,1) = std::complex<double>(pars(271), pars(272));
   UM(2,2) = std::complex<double>(pars(273), pars(274));
   UP(0,0) = std::complex<double>(pars(275), pars(276));
   UP(0,1) = std::complex<double>(pars(277), pars(278));
   UP(0,2) = std::complex<double>(pars(279), pars(280));
   UP(1,0) = std::complex<double>(pars(281), pars(282));
   UP(1,1) = std::complex<double>(pars(283), pars(284));
   UP(1,2) = std::complex<double>(pars(285), pars(286));
   UP(2,0) = std::complex<double>(pars(287), pars(288));
   UP(2,1) = std::complex<double>(pars(289), pars(290));
   UP(2,2) = std::complex<double>(pars(291), pars(292));
   ZEL(0,0) = std::complex<double>(pars(293), pars(294));
   ZEL(0,1) = std::complex<double>(pars(295), pars(296));
   ZEL(0,2) = std::complex<double>(pars(297), pars(298));
   ZEL(1,0) = std::complex<double>(pars(299), pars(300));
   ZEL(1,1) = std::complex<double>(pars(301), pars(302));
   ZEL(1,2) = std::complex<double>(pars(303), pars(304));
   ZEL(2,0) = std::complex<double>(pars(305), pars(306));
   ZEL(2,1) = std::complex<double>(pars(307), pars(308));
   ZEL(2,2) = std::complex<double>(pars(309), pars(310));
   ZER(0,0) = std::complex<double>(pars(311), pars(312));
   ZER(0,1) = std::complex<double>(pars(313), pars(314));
   ZER(0,2) = std::complex<double>(pars(315), pars(316));
   ZER(1,0) = std::complex<double>(pars(317), pars(318));
   ZER(1,1) = std::complex<double>(pars(319), pars(320));
   ZER(1,2) = std::complex<double>(pars(321), pars(322));
   ZER(2,0) = std::complex<double>(pars(323), pars(324));
   ZER(2,1) = std::complex<double>(pars(325), pars(326));
   ZER(2,2) = std::complex<double>(pars(327), pars(328));
   ZDL(0,0) = std::complex<double>(pars(329), pars(330));
   ZDL(0,1) = std::complex<double>(pars(331), pars(332));
   ZDL(0,2) = std::complex<double>(pars(333), pars(334));
   ZDL(1,0) = std::complex<double>(pars(335), pars(336));
   ZDL(1,1) = std::complex<double>(pars(337), pars(338));
   ZDL(1,2) = std::complex<double>(pars(339), pars(340));
   ZDL(2,0) = std::complex<double>(pars(341), pars(342));
   ZDL(2,1) = std::complex<double>(pars(343), pars(344));
   ZDL(2,2) = std::complex<double>(pars(345), pars(346));
   ZDR(0,0) = std::complex<double>(pars(347), pars(348));
   ZDR(0,1) = std::complex<double>(pars(349), pars(350));
   ZDR(0,2) = std::complex<double>(pars(351), pars(352));
   ZDR(1,0) = std::complex<double>(pars(353), pars(354));
   ZDR(1,1) = std::complex<double>(pars(355), pars(356));
   ZDR(1,2) = std::complex<double>(pars(357), pars(358));
   ZDR(2,0) = std::complex<double>(pars(359), pars(360));
   ZDR(2,1) = std::complex<double>(pars(361), pars(362));
   ZDR(2,2) = std::complex<double>(pars(363), pars(364));
   ZUL(0,0) = std::complex<double>(pars(365), pars(366));
   ZUL(0,1) = std::complex<double>(pars(367), pars(368));
   ZUL(0,2) = std::complex<double>(pars(369), pars(370));
   ZUL(1,0) = std::complex<double>(pars(371), pars(372));
   ZUL(1,1) = std::complex<double>(pars(373), pars(374));
   ZUL(1,2) = std::complex<double>(pars(375), pars(376));
   ZUL(2,0) = std::complex<double>(pars(377), pars(378));
   ZUL(2,1) = std::complex<double>(pars(379), pars(380));
   ZUL(2,2) = std::complex<double>(pars(381), pars(382));
   ZUR(0,0) = std::complex<double>(pars(383), pars(384));
   ZUR(0,1) = std::complex<double>(pars(385), pars(386));
   ZUR(0,2) = std::complex<double>(pars(387), pars(388));
   ZUR(1,0) = std::complex<double>(pars(389), pars(390));
   ZUR(1,1) = std::complex<double>(pars(391), pars(392));
   ZUR(1,2) = std::complex<double>(pars(393), pars(394));
   ZUR(2,0) = std::complex<double>(pars(395), pars(396));
   ZUR(2,1) = std::complex<double>(pars(397), pars(398));
   ZUR(2,2) = std::complex<double>(pars(399), pars(400));
   ZZ(0,0) = pars(401);
   ZZ(0,1) = pars(402);
   ZZ(1,0) = pars(403);
   ZZ(1,1) = pars(404);

}

Eigen::ArrayXd TMSSM_physical::get_masses() const
{
   Eigen::ArrayXd pars(56);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MSd(0);
   pars(6) = MSd(1);
   pars(7) = MSd(2);
   pars(8) = MSd(3);
   pars(9) = MSd(4);
   pars(10) = MSd(5);
   pars(11) = MSv(0);
   pars(12) = MSv(1);
   pars(13) = MSv(2);
   pars(14) = MSu(0);
   pars(15) = MSu(1);
   pars(16) = MSu(2);
   pars(17) = MSu(3);
   pars(18) = MSu(4);
   pars(19) = MSu(5);
   pars(20) = MSe(0);
   pars(21) = MSe(1);
   pars(22) = MSe(2);
   pars(23) = MSe(3);
   pars(24) = MSe(4);
   pars(25) = MSe(5);
   pars(26) = Mhh(0);
   pars(27) = Mhh(1);
   pars(28) = Mhh(2);
   pars(29) = MAh(0);
   pars(30) = MAh(1);
   pars(31) = MAh(2);
   pars(32) = MHpm(0);
   pars(33) = MHpm(1);
   pars(34) = MHpm(2);
   pars(35) = MHpm(3);
   pars(36) = MChi(0);
   pars(37) = MChi(1);
   pars(38) = MChi(2);
   pars(39) = MChi(3);
   pars(40) = MChi(4);
   pars(41) = MCha(0);
   pars(42) = MCha(1);
   pars(43) = MCha(2);
   pars(44) = MFe(0);
   pars(45) = MFe(1);
   pars(46) = MFe(2);
   pars(47) = MFd(0);
   pars(48) = MFd(1);
   pars(49) = MFd(2);
   pars(50) = MFu(0);
   pars(51) = MFu(1);
   pars(52) = MFu(2);
   pars(53) = MVWm;
   pars(54) = MVP;
   pars(55) = MVZ;

   return pars;
}

void TMSSM_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MFv(0) = pars(2);
   MFv(1) = pars(3);
   MFv(2) = pars(4);
   MSd(0) = pars(5);
   MSd(1) = pars(6);
   MSd(2) = pars(7);
   MSd(3) = pars(8);
   MSd(4) = pars(9);
   MSd(5) = pars(10);
   MSv(0) = pars(11);
   MSv(1) = pars(12);
   MSv(2) = pars(13);
   MSu(0) = pars(14);
   MSu(1) = pars(15);
   MSu(2) = pars(16);
   MSu(3) = pars(17);
   MSu(4) = pars(18);
   MSu(5) = pars(19);
   MSe(0) = pars(20);
   MSe(1) = pars(21);
   MSe(2) = pars(22);
   MSe(3) = pars(23);
   MSe(4) = pars(24);
   MSe(5) = pars(25);
   Mhh(0) = pars(26);
   Mhh(1) = pars(27);
   Mhh(2) = pars(28);
   MAh(0) = pars(29);
   MAh(1) = pars(30);
   MAh(2) = pars(31);
   MHpm(0) = pars(32);
   MHpm(1) = pars(33);
   MHpm(2) = pars(34);
   MHpm(3) = pars(35);
   MChi(0) = pars(36);
   MChi(1) = pars(37);
   MChi(2) = pars(38);
   MChi(3) = pars(39);
   MChi(4) = pars(40);
   MCha(0) = pars(41);
   MCha(1) = pars(42);
   MCha(2) = pars(43);
   MFe(0) = pars(44);
   MFe(1) = pars(45);
   MFe(2) = pars(46);
   MFd(0) = pars(47);
   MFd(1) = pars(48);
   MFd(2) = pars(49);
   MFu(0) = pars(50);
   MFu(1) = pars(51);
   MFu(2) = pars(52);
   MVWm = pars(53);
   MVP = pars(54);
   MVZ = pars(55);

}

void TMSSM_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha = " << MCha.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MVWm = " << MVWm << '\n';

   ostr << "----------------------------------------\n"
           "pole mass mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "ZD = " << ZD << '\n';
   ostr << "ZV = " << ZV << '\n';
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
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
   ostr << "ZZ = " << ZZ << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const TMSSM_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
