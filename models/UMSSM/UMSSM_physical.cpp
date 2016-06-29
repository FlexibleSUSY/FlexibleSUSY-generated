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

// File generated at Wed 29 Jun 2016 12:21:23

#include "UMSSM_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

UMSSM_physical::UMSSM_physical()
   :
    MVG(0), MGlu(0), MSd(Eigen::Array<double,6,1>::Zero()), MSv(Eigen::Array<
       double,6,1>::Zero()), MSu(Eigen::Array<double,6,1>::Zero()), MSe(
       Eigen::Array<double,6,1>::Zero()), Mhh(Eigen::Array<double,3,1>::Zero()),
       MAh(Eigen::Array<double,3,1>::Zero()), MHpm(Eigen::Array<double,2,1>::Zero(
       )), MChi(Eigen::Array<double,6,1>::Zero()), MFv(Eigen::Array<double,3,1>
       ::Zero()), MCha(Eigen::Array<double,2,1>::Zero()), MFe(Eigen::Array<double,
       3,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<
       double,3,1>::Zero()), MVWm(0), MVP(0), MVZ(0), MVZp(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,6,6>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,3,3>::Zero()), ZA(Eigen::Matrix<double,3,
      3>::Zero()), ZP(Eigen::Matrix<double,2,2>::Zero()), ZN(Eigen::Matrix<
      std::complex<double>,6,6>::Zero()), ZVL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZVR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), UM(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP(Eigen::Matrix<
      std::complex<double>,2,2>::Zero()), ZEL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZER(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDL(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDR(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZUL(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZUR(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZZ(
      Eigen::Matrix<double,3,3>::Zero())

{
}

void UMSSM_physical::clear()
{
   MVG = 0.;
   MGlu = 0.;
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,6,1>::Zero();
   ZV = Eigen::Matrix<double,6,6>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,3,1>::Zero();
   ZH = Eigen::Matrix<double,3,3>::Zero();
   MAh = Eigen::Matrix<double,3,1>::Zero();
   ZA = Eigen::Matrix<double,3,3>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,6,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,6,6>::Zero();
   MFv = Eigen::Matrix<double,3,1>::Zero();
   ZVL = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   ZVR = Eigen::Matrix<std::complex<double>,3,3>::Zero();
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
void UMSSM_physical::convert_to_hk()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void UMSSM_physical::convert_to_slha()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));

}

Eigen::ArrayXd UMSSM_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(465);

   pars(58) = ZD(0,0);
   pars(59) = ZD(0,1);
   pars(60) = ZD(0,2);
   pars(61) = ZD(0,3);
   pars(62) = ZD(0,4);
   pars(63) = ZD(0,5);
   pars(64) = ZD(1,0);
   pars(65) = ZD(1,1);
   pars(66) = ZD(1,2);
   pars(67) = ZD(1,3);
   pars(68) = ZD(1,4);
   pars(69) = ZD(1,5);
   pars(70) = ZD(2,0);
   pars(71) = ZD(2,1);
   pars(72) = ZD(2,2);
   pars(73) = ZD(2,3);
   pars(74) = ZD(2,4);
   pars(75) = ZD(2,5);
   pars(76) = ZD(3,0);
   pars(77) = ZD(3,1);
   pars(78) = ZD(3,2);
   pars(79) = ZD(3,3);
   pars(80) = ZD(3,4);
   pars(81) = ZD(3,5);
   pars(82) = ZD(4,0);
   pars(83) = ZD(4,1);
   pars(84) = ZD(4,2);
   pars(85) = ZD(4,3);
   pars(86) = ZD(4,4);
   pars(87) = ZD(4,5);
   pars(88) = ZD(5,0);
   pars(89) = ZD(5,1);
   pars(90) = ZD(5,2);
   pars(91) = ZD(5,3);
   pars(92) = ZD(5,4);
   pars(93) = ZD(5,5);
   pars(94) = ZV(0,0);
   pars(95) = ZV(0,1);
   pars(96) = ZV(0,2);
   pars(97) = ZV(0,3);
   pars(98) = ZV(0,4);
   pars(99) = ZV(0,5);
   pars(100) = ZV(1,0);
   pars(101) = ZV(1,1);
   pars(102) = ZV(1,2);
   pars(103) = ZV(1,3);
   pars(104) = ZV(1,4);
   pars(105) = ZV(1,5);
   pars(106) = ZV(2,0);
   pars(107) = ZV(2,1);
   pars(108) = ZV(2,2);
   pars(109) = ZV(2,3);
   pars(110) = ZV(2,4);
   pars(111) = ZV(2,5);
   pars(112) = ZV(3,0);
   pars(113) = ZV(3,1);
   pars(114) = ZV(3,2);
   pars(115) = ZV(3,3);
   pars(116) = ZV(3,4);
   pars(117) = ZV(3,5);
   pars(118) = ZV(4,0);
   pars(119) = ZV(4,1);
   pars(120) = ZV(4,2);
   pars(121) = ZV(4,3);
   pars(122) = ZV(4,4);
   pars(123) = ZV(4,5);
   pars(124) = ZV(5,0);
   pars(125) = ZV(5,1);
   pars(126) = ZV(5,2);
   pars(127) = ZV(5,3);
   pars(128) = ZV(5,4);
   pars(129) = ZV(5,5);
   pars(130) = ZU(0,0);
   pars(131) = ZU(0,1);
   pars(132) = ZU(0,2);
   pars(133) = ZU(0,3);
   pars(134) = ZU(0,4);
   pars(135) = ZU(0,5);
   pars(136) = ZU(1,0);
   pars(137) = ZU(1,1);
   pars(138) = ZU(1,2);
   pars(139) = ZU(1,3);
   pars(140) = ZU(1,4);
   pars(141) = ZU(1,5);
   pars(142) = ZU(2,0);
   pars(143) = ZU(2,1);
   pars(144) = ZU(2,2);
   pars(145) = ZU(2,3);
   pars(146) = ZU(2,4);
   pars(147) = ZU(2,5);
   pars(148) = ZU(3,0);
   pars(149) = ZU(3,1);
   pars(150) = ZU(3,2);
   pars(151) = ZU(3,3);
   pars(152) = ZU(3,4);
   pars(153) = ZU(3,5);
   pars(154) = ZU(4,0);
   pars(155) = ZU(4,1);
   pars(156) = ZU(4,2);
   pars(157) = ZU(4,3);
   pars(158) = ZU(4,4);
   pars(159) = ZU(4,5);
   pars(160) = ZU(5,0);
   pars(161) = ZU(5,1);
   pars(162) = ZU(5,2);
   pars(163) = ZU(5,3);
   pars(164) = ZU(5,4);
   pars(165) = ZU(5,5);
   pars(166) = ZE(0,0);
   pars(167) = ZE(0,1);
   pars(168) = ZE(0,2);
   pars(169) = ZE(0,3);
   pars(170) = ZE(0,4);
   pars(171) = ZE(0,5);
   pars(172) = ZE(1,0);
   pars(173) = ZE(1,1);
   pars(174) = ZE(1,2);
   pars(175) = ZE(1,3);
   pars(176) = ZE(1,4);
   pars(177) = ZE(1,5);
   pars(178) = ZE(2,0);
   pars(179) = ZE(2,1);
   pars(180) = ZE(2,2);
   pars(181) = ZE(2,3);
   pars(182) = ZE(2,4);
   pars(183) = ZE(2,5);
   pars(184) = ZE(3,0);
   pars(185) = ZE(3,1);
   pars(186) = ZE(3,2);
   pars(187) = ZE(3,3);
   pars(188) = ZE(3,4);
   pars(189) = ZE(3,5);
   pars(190) = ZE(4,0);
   pars(191) = ZE(4,1);
   pars(192) = ZE(4,2);
   pars(193) = ZE(4,3);
   pars(194) = ZE(4,4);
   pars(195) = ZE(4,5);
   pars(196) = ZE(5,0);
   pars(197) = ZE(5,1);
   pars(198) = ZE(5,2);
   pars(199) = ZE(5,3);
   pars(200) = ZE(5,4);
   pars(201) = ZE(5,5);
   pars(202) = ZH(0,0);
   pars(203) = ZH(0,1);
   pars(204) = ZH(0,2);
   pars(205) = ZH(1,0);
   pars(206) = ZH(1,1);
   pars(207) = ZH(1,2);
   pars(208) = ZH(2,0);
   pars(209) = ZH(2,1);
   pars(210) = ZH(2,2);
   pars(211) = ZA(0,0);
   pars(212) = ZA(0,1);
   pars(213) = ZA(0,2);
   pars(214) = ZA(1,0);
   pars(215) = ZA(1,1);
   pars(216) = ZA(1,2);
   pars(217) = ZA(2,0);
   pars(218) = ZA(2,1);
   pars(219) = ZA(2,2);
   pars(220) = ZP(0,0);
   pars(221) = ZP(0,1);
   pars(222) = ZP(1,0);
   pars(223) = ZP(1,1);
   pars(224) = Re(ZN(0,0));
   pars(225) = Im(ZN(0,0));
   pars(226) = Re(ZN(0,1));
   pars(227) = Im(ZN(0,1));
   pars(228) = Re(ZN(0,2));
   pars(229) = Im(ZN(0,2));
   pars(230) = Re(ZN(0,3));
   pars(231) = Im(ZN(0,3));
   pars(232) = Re(ZN(0,4));
   pars(233) = Im(ZN(0,4));
   pars(234) = Re(ZN(0,5));
   pars(235) = Im(ZN(0,5));
   pars(236) = Re(ZN(1,0));
   pars(237) = Im(ZN(1,0));
   pars(238) = Re(ZN(1,1));
   pars(239) = Im(ZN(1,1));
   pars(240) = Re(ZN(1,2));
   pars(241) = Im(ZN(1,2));
   pars(242) = Re(ZN(1,3));
   pars(243) = Im(ZN(1,3));
   pars(244) = Re(ZN(1,4));
   pars(245) = Im(ZN(1,4));
   pars(246) = Re(ZN(1,5));
   pars(247) = Im(ZN(1,5));
   pars(248) = Re(ZN(2,0));
   pars(249) = Im(ZN(2,0));
   pars(250) = Re(ZN(2,1));
   pars(251) = Im(ZN(2,1));
   pars(252) = Re(ZN(2,2));
   pars(253) = Im(ZN(2,2));
   pars(254) = Re(ZN(2,3));
   pars(255) = Im(ZN(2,3));
   pars(256) = Re(ZN(2,4));
   pars(257) = Im(ZN(2,4));
   pars(258) = Re(ZN(2,5));
   pars(259) = Im(ZN(2,5));
   pars(260) = Re(ZN(3,0));
   pars(261) = Im(ZN(3,0));
   pars(262) = Re(ZN(3,1));
   pars(263) = Im(ZN(3,1));
   pars(264) = Re(ZN(3,2));
   pars(265) = Im(ZN(3,2));
   pars(266) = Re(ZN(3,3));
   pars(267) = Im(ZN(3,3));
   pars(268) = Re(ZN(3,4));
   pars(269) = Im(ZN(3,4));
   pars(270) = Re(ZN(3,5));
   pars(271) = Im(ZN(3,5));
   pars(272) = Re(ZN(4,0));
   pars(273) = Im(ZN(4,0));
   pars(274) = Re(ZN(4,1));
   pars(275) = Im(ZN(4,1));
   pars(276) = Re(ZN(4,2));
   pars(277) = Im(ZN(4,2));
   pars(278) = Re(ZN(4,3));
   pars(279) = Im(ZN(4,3));
   pars(280) = Re(ZN(4,4));
   pars(281) = Im(ZN(4,4));
   pars(282) = Re(ZN(4,5));
   pars(283) = Im(ZN(4,5));
   pars(284) = Re(ZN(5,0));
   pars(285) = Im(ZN(5,0));
   pars(286) = Re(ZN(5,1));
   pars(287) = Im(ZN(5,1));
   pars(288) = Re(ZN(5,2));
   pars(289) = Im(ZN(5,2));
   pars(290) = Re(ZN(5,3));
   pars(291) = Im(ZN(5,3));
   pars(292) = Re(ZN(5,4));
   pars(293) = Im(ZN(5,4));
   pars(294) = Re(ZN(5,5));
   pars(295) = Im(ZN(5,5));
   pars(296) = Re(ZVL(0,0));
   pars(297) = Im(ZVL(0,0));
   pars(298) = Re(ZVL(0,1));
   pars(299) = Im(ZVL(0,1));
   pars(300) = Re(ZVL(0,2));
   pars(301) = Im(ZVL(0,2));
   pars(302) = Re(ZVL(1,0));
   pars(303) = Im(ZVL(1,0));
   pars(304) = Re(ZVL(1,1));
   pars(305) = Im(ZVL(1,1));
   pars(306) = Re(ZVL(1,2));
   pars(307) = Im(ZVL(1,2));
   pars(308) = Re(ZVL(2,0));
   pars(309) = Im(ZVL(2,0));
   pars(310) = Re(ZVL(2,1));
   pars(311) = Im(ZVL(2,1));
   pars(312) = Re(ZVL(2,2));
   pars(313) = Im(ZVL(2,2));
   pars(314) = Re(ZVR(0,0));
   pars(315) = Im(ZVR(0,0));
   pars(316) = Re(ZVR(0,1));
   pars(317) = Im(ZVR(0,1));
   pars(318) = Re(ZVR(0,2));
   pars(319) = Im(ZVR(0,2));
   pars(320) = Re(ZVR(1,0));
   pars(321) = Im(ZVR(1,0));
   pars(322) = Re(ZVR(1,1));
   pars(323) = Im(ZVR(1,1));
   pars(324) = Re(ZVR(1,2));
   pars(325) = Im(ZVR(1,2));
   pars(326) = Re(ZVR(2,0));
   pars(327) = Im(ZVR(2,0));
   pars(328) = Re(ZVR(2,1));
   pars(329) = Im(ZVR(2,1));
   pars(330) = Re(ZVR(2,2));
   pars(331) = Im(ZVR(2,2));
   pars(332) = Re(UM(0,0));
   pars(333) = Im(UM(0,0));
   pars(334) = Re(UM(0,1));
   pars(335) = Im(UM(0,1));
   pars(336) = Re(UM(1,0));
   pars(337) = Im(UM(1,0));
   pars(338) = Re(UM(1,1));
   pars(339) = Im(UM(1,1));
   pars(340) = Re(UP(0,0));
   pars(341) = Im(UP(0,0));
   pars(342) = Re(UP(0,1));
   pars(343) = Im(UP(0,1));
   pars(344) = Re(UP(1,0));
   pars(345) = Im(UP(1,0));
   pars(346) = Re(UP(1,1));
   pars(347) = Im(UP(1,1));
   pars(348) = Re(ZEL(0,0));
   pars(349) = Im(ZEL(0,0));
   pars(350) = Re(ZEL(0,1));
   pars(351) = Im(ZEL(0,1));
   pars(352) = Re(ZEL(0,2));
   pars(353) = Im(ZEL(0,2));
   pars(354) = Re(ZEL(1,0));
   pars(355) = Im(ZEL(1,0));
   pars(356) = Re(ZEL(1,1));
   pars(357) = Im(ZEL(1,1));
   pars(358) = Re(ZEL(1,2));
   pars(359) = Im(ZEL(1,2));
   pars(360) = Re(ZEL(2,0));
   pars(361) = Im(ZEL(2,0));
   pars(362) = Re(ZEL(2,1));
   pars(363) = Im(ZEL(2,1));
   pars(364) = Re(ZEL(2,2));
   pars(365) = Im(ZEL(2,2));
   pars(366) = Re(ZER(0,0));
   pars(367) = Im(ZER(0,0));
   pars(368) = Re(ZER(0,1));
   pars(369) = Im(ZER(0,1));
   pars(370) = Re(ZER(0,2));
   pars(371) = Im(ZER(0,2));
   pars(372) = Re(ZER(1,0));
   pars(373) = Im(ZER(1,0));
   pars(374) = Re(ZER(1,1));
   pars(375) = Im(ZER(1,1));
   pars(376) = Re(ZER(1,2));
   pars(377) = Im(ZER(1,2));
   pars(378) = Re(ZER(2,0));
   pars(379) = Im(ZER(2,0));
   pars(380) = Re(ZER(2,1));
   pars(381) = Im(ZER(2,1));
   pars(382) = Re(ZER(2,2));
   pars(383) = Im(ZER(2,2));
   pars(384) = Re(ZDL(0,0));
   pars(385) = Im(ZDL(0,0));
   pars(386) = Re(ZDL(0,1));
   pars(387) = Im(ZDL(0,1));
   pars(388) = Re(ZDL(0,2));
   pars(389) = Im(ZDL(0,2));
   pars(390) = Re(ZDL(1,0));
   pars(391) = Im(ZDL(1,0));
   pars(392) = Re(ZDL(1,1));
   pars(393) = Im(ZDL(1,1));
   pars(394) = Re(ZDL(1,2));
   pars(395) = Im(ZDL(1,2));
   pars(396) = Re(ZDL(2,0));
   pars(397) = Im(ZDL(2,0));
   pars(398) = Re(ZDL(2,1));
   pars(399) = Im(ZDL(2,1));
   pars(400) = Re(ZDL(2,2));
   pars(401) = Im(ZDL(2,2));
   pars(402) = Re(ZDR(0,0));
   pars(403) = Im(ZDR(0,0));
   pars(404) = Re(ZDR(0,1));
   pars(405) = Im(ZDR(0,1));
   pars(406) = Re(ZDR(0,2));
   pars(407) = Im(ZDR(0,2));
   pars(408) = Re(ZDR(1,0));
   pars(409) = Im(ZDR(1,0));
   pars(410) = Re(ZDR(1,1));
   pars(411) = Im(ZDR(1,1));
   pars(412) = Re(ZDR(1,2));
   pars(413) = Im(ZDR(1,2));
   pars(414) = Re(ZDR(2,0));
   pars(415) = Im(ZDR(2,0));
   pars(416) = Re(ZDR(2,1));
   pars(417) = Im(ZDR(2,1));
   pars(418) = Re(ZDR(2,2));
   pars(419) = Im(ZDR(2,2));
   pars(420) = Re(ZUL(0,0));
   pars(421) = Im(ZUL(0,0));
   pars(422) = Re(ZUL(0,1));
   pars(423) = Im(ZUL(0,1));
   pars(424) = Re(ZUL(0,2));
   pars(425) = Im(ZUL(0,2));
   pars(426) = Re(ZUL(1,0));
   pars(427) = Im(ZUL(1,0));
   pars(428) = Re(ZUL(1,1));
   pars(429) = Im(ZUL(1,1));
   pars(430) = Re(ZUL(1,2));
   pars(431) = Im(ZUL(1,2));
   pars(432) = Re(ZUL(2,0));
   pars(433) = Im(ZUL(2,0));
   pars(434) = Re(ZUL(2,1));
   pars(435) = Im(ZUL(2,1));
   pars(436) = Re(ZUL(2,2));
   pars(437) = Im(ZUL(2,2));
   pars(438) = Re(ZUR(0,0));
   pars(439) = Im(ZUR(0,0));
   pars(440) = Re(ZUR(0,1));
   pars(441) = Im(ZUR(0,1));
   pars(442) = Re(ZUR(0,2));
   pars(443) = Im(ZUR(0,2));
   pars(444) = Re(ZUR(1,0));
   pars(445) = Im(ZUR(1,0));
   pars(446) = Re(ZUR(1,1));
   pars(447) = Im(ZUR(1,1));
   pars(448) = Re(ZUR(1,2));
   pars(449) = Im(ZUR(1,2));
   pars(450) = Re(ZUR(2,0));
   pars(451) = Im(ZUR(2,0));
   pars(452) = Re(ZUR(2,1));
   pars(453) = Im(ZUR(2,1));
   pars(454) = Re(ZUR(2,2));
   pars(455) = Im(ZUR(2,2));
   pars(456) = ZZ(0,0);
   pars(457) = ZZ(0,1);
   pars(458) = ZZ(0,2);
   pars(459) = ZZ(1,0);
   pars(460) = ZZ(1,1);
   pars(461) = ZZ(1,2);
   pars(462) = ZZ(2,0);
   pars(463) = ZZ(2,1);
   pars(464) = ZZ(2,2);


   return pars;
}

void UMSSM_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

   ZD(0,0) = pars(58);
   ZD(0,1) = pars(59);
   ZD(0,2) = pars(60);
   ZD(0,3) = pars(61);
   ZD(0,4) = pars(62);
   ZD(0,5) = pars(63);
   ZD(1,0) = pars(64);
   ZD(1,1) = pars(65);
   ZD(1,2) = pars(66);
   ZD(1,3) = pars(67);
   ZD(1,4) = pars(68);
   ZD(1,5) = pars(69);
   ZD(2,0) = pars(70);
   ZD(2,1) = pars(71);
   ZD(2,2) = pars(72);
   ZD(2,3) = pars(73);
   ZD(2,4) = pars(74);
   ZD(2,5) = pars(75);
   ZD(3,0) = pars(76);
   ZD(3,1) = pars(77);
   ZD(3,2) = pars(78);
   ZD(3,3) = pars(79);
   ZD(3,4) = pars(80);
   ZD(3,5) = pars(81);
   ZD(4,0) = pars(82);
   ZD(4,1) = pars(83);
   ZD(4,2) = pars(84);
   ZD(4,3) = pars(85);
   ZD(4,4) = pars(86);
   ZD(4,5) = pars(87);
   ZD(5,0) = pars(88);
   ZD(5,1) = pars(89);
   ZD(5,2) = pars(90);
   ZD(5,3) = pars(91);
   ZD(5,4) = pars(92);
   ZD(5,5) = pars(93);
   ZV(0,0) = pars(94);
   ZV(0,1) = pars(95);
   ZV(0,2) = pars(96);
   ZV(0,3) = pars(97);
   ZV(0,4) = pars(98);
   ZV(0,5) = pars(99);
   ZV(1,0) = pars(100);
   ZV(1,1) = pars(101);
   ZV(1,2) = pars(102);
   ZV(1,3) = pars(103);
   ZV(1,4) = pars(104);
   ZV(1,5) = pars(105);
   ZV(2,0) = pars(106);
   ZV(2,1) = pars(107);
   ZV(2,2) = pars(108);
   ZV(2,3) = pars(109);
   ZV(2,4) = pars(110);
   ZV(2,5) = pars(111);
   ZV(3,0) = pars(112);
   ZV(3,1) = pars(113);
   ZV(3,2) = pars(114);
   ZV(3,3) = pars(115);
   ZV(3,4) = pars(116);
   ZV(3,5) = pars(117);
   ZV(4,0) = pars(118);
   ZV(4,1) = pars(119);
   ZV(4,2) = pars(120);
   ZV(4,3) = pars(121);
   ZV(4,4) = pars(122);
   ZV(4,5) = pars(123);
   ZV(5,0) = pars(124);
   ZV(5,1) = pars(125);
   ZV(5,2) = pars(126);
   ZV(5,3) = pars(127);
   ZV(5,4) = pars(128);
   ZV(5,5) = pars(129);
   ZU(0,0) = pars(130);
   ZU(0,1) = pars(131);
   ZU(0,2) = pars(132);
   ZU(0,3) = pars(133);
   ZU(0,4) = pars(134);
   ZU(0,5) = pars(135);
   ZU(1,0) = pars(136);
   ZU(1,1) = pars(137);
   ZU(1,2) = pars(138);
   ZU(1,3) = pars(139);
   ZU(1,4) = pars(140);
   ZU(1,5) = pars(141);
   ZU(2,0) = pars(142);
   ZU(2,1) = pars(143);
   ZU(2,2) = pars(144);
   ZU(2,3) = pars(145);
   ZU(2,4) = pars(146);
   ZU(2,5) = pars(147);
   ZU(3,0) = pars(148);
   ZU(3,1) = pars(149);
   ZU(3,2) = pars(150);
   ZU(3,3) = pars(151);
   ZU(3,4) = pars(152);
   ZU(3,5) = pars(153);
   ZU(4,0) = pars(154);
   ZU(4,1) = pars(155);
   ZU(4,2) = pars(156);
   ZU(4,3) = pars(157);
   ZU(4,4) = pars(158);
   ZU(4,5) = pars(159);
   ZU(5,0) = pars(160);
   ZU(5,1) = pars(161);
   ZU(5,2) = pars(162);
   ZU(5,3) = pars(163);
   ZU(5,4) = pars(164);
   ZU(5,5) = pars(165);
   ZE(0,0) = pars(166);
   ZE(0,1) = pars(167);
   ZE(0,2) = pars(168);
   ZE(0,3) = pars(169);
   ZE(0,4) = pars(170);
   ZE(0,5) = pars(171);
   ZE(1,0) = pars(172);
   ZE(1,1) = pars(173);
   ZE(1,2) = pars(174);
   ZE(1,3) = pars(175);
   ZE(1,4) = pars(176);
   ZE(1,5) = pars(177);
   ZE(2,0) = pars(178);
   ZE(2,1) = pars(179);
   ZE(2,2) = pars(180);
   ZE(2,3) = pars(181);
   ZE(2,4) = pars(182);
   ZE(2,5) = pars(183);
   ZE(3,0) = pars(184);
   ZE(3,1) = pars(185);
   ZE(3,2) = pars(186);
   ZE(3,3) = pars(187);
   ZE(3,4) = pars(188);
   ZE(3,5) = pars(189);
   ZE(4,0) = pars(190);
   ZE(4,1) = pars(191);
   ZE(4,2) = pars(192);
   ZE(4,3) = pars(193);
   ZE(4,4) = pars(194);
   ZE(4,5) = pars(195);
   ZE(5,0) = pars(196);
   ZE(5,1) = pars(197);
   ZE(5,2) = pars(198);
   ZE(5,3) = pars(199);
   ZE(5,4) = pars(200);
   ZE(5,5) = pars(201);
   ZH(0,0) = pars(202);
   ZH(0,1) = pars(203);
   ZH(0,2) = pars(204);
   ZH(1,0) = pars(205);
   ZH(1,1) = pars(206);
   ZH(1,2) = pars(207);
   ZH(2,0) = pars(208);
   ZH(2,1) = pars(209);
   ZH(2,2) = pars(210);
   ZA(0,0) = pars(211);
   ZA(0,1) = pars(212);
   ZA(0,2) = pars(213);
   ZA(1,0) = pars(214);
   ZA(1,1) = pars(215);
   ZA(1,2) = pars(216);
   ZA(2,0) = pars(217);
   ZA(2,1) = pars(218);
   ZA(2,2) = pars(219);
   ZP(0,0) = pars(220);
   ZP(0,1) = pars(221);
   ZP(1,0) = pars(222);
   ZP(1,1) = pars(223);
   ZN(0,0) = std::complex<double>(pars(224), pars(225));
   ZN(0,1) = std::complex<double>(pars(226), pars(227));
   ZN(0,2) = std::complex<double>(pars(228), pars(229));
   ZN(0,3) = std::complex<double>(pars(230), pars(231));
   ZN(0,4) = std::complex<double>(pars(232), pars(233));
   ZN(0,5) = std::complex<double>(pars(234), pars(235));
   ZN(1,0) = std::complex<double>(pars(236), pars(237));
   ZN(1,1) = std::complex<double>(pars(238), pars(239));
   ZN(1,2) = std::complex<double>(pars(240), pars(241));
   ZN(1,3) = std::complex<double>(pars(242), pars(243));
   ZN(1,4) = std::complex<double>(pars(244), pars(245));
   ZN(1,5) = std::complex<double>(pars(246), pars(247));
   ZN(2,0) = std::complex<double>(pars(248), pars(249));
   ZN(2,1) = std::complex<double>(pars(250), pars(251));
   ZN(2,2) = std::complex<double>(pars(252), pars(253));
   ZN(2,3) = std::complex<double>(pars(254), pars(255));
   ZN(2,4) = std::complex<double>(pars(256), pars(257));
   ZN(2,5) = std::complex<double>(pars(258), pars(259));
   ZN(3,0) = std::complex<double>(pars(260), pars(261));
   ZN(3,1) = std::complex<double>(pars(262), pars(263));
   ZN(3,2) = std::complex<double>(pars(264), pars(265));
   ZN(3,3) = std::complex<double>(pars(266), pars(267));
   ZN(3,4) = std::complex<double>(pars(268), pars(269));
   ZN(3,5) = std::complex<double>(pars(270), pars(271));
   ZN(4,0) = std::complex<double>(pars(272), pars(273));
   ZN(4,1) = std::complex<double>(pars(274), pars(275));
   ZN(4,2) = std::complex<double>(pars(276), pars(277));
   ZN(4,3) = std::complex<double>(pars(278), pars(279));
   ZN(4,4) = std::complex<double>(pars(280), pars(281));
   ZN(4,5) = std::complex<double>(pars(282), pars(283));
   ZN(5,0) = std::complex<double>(pars(284), pars(285));
   ZN(5,1) = std::complex<double>(pars(286), pars(287));
   ZN(5,2) = std::complex<double>(pars(288), pars(289));
   ZN(5,3) = std::complex<double>(pars(290), pars(291));
   ZN(5,4) = std::complex<double>(pars(292), pars(293));
   ZN(5,5) = std::complex<double>(pars(294), pars(295));
   ZVL(0,0) = std::complex<double>(pars(296), pars(297));
   ZVL(0,1) = std::complex<double>(pars(298), pars(299));
   ZVL(0,2) = std::complex<double>(pars(300), pars(301));
   ZVL(1,0) = std::complex<double>(pars(302), pars(303));
   ZVL(1,1) = std::complex<double>(pars(304), pars(305));
   ZVL(1,2) = std::complex<double>(pars(306), pars(307));
   ZVL(2,0) = std::complex<double>(pars(308), pars(309));
   ZVL(2,1) = std::complex<double>(pars(310), pars(311));
   ZVL(2,2) = std::complex<double>(pars(312), pars(313));
   ZVR(0,0) = std::complex<double>(pars(314), pars(315));
   ZVR(0,1) = std::complex<double>(pars(316), pars(317));
   ZVR(0,2) = std::complex<double>(pars(318), pars(319));
   ZVR(1,0) = std::complex<double>(pars(320), pars(321));
   ZVR(1,1) = std::complex<double>(pars(322), pars(323));
   ZVR(1,2) = std::complex<double>(pars(324), pars(325));
   ZVR(2,0) = std::complex<double>(pars(326), pars(327));
   ZVR(2,1) = std::complex<double>(pars(328), pars(329));
   ZVR(2,2) = std::complex<double>(pars(330), pars(331));
   UM(0,0) = std::complex<double>(pars(332), pars(333));
   UM(0,1) = std::complex<double>(pars(334), pars(335));
   UM(1,0) = std::complex<double>(pars(336), pars(337));
   UM(1,1) = std::complex<double>(pars(338), pars(339));
   UP(0,0) = std::complex<double>(pars(340), pars(341));
   UP(0,1) = std::complex<double>(pars(342), pars(343));
   UP(1,0) = std::complex<double>(pars(344), pars(345));
   UP(1,1) = std::complex<double>(pars(346), pars(347));
   ZEL(0,0) = std::complex<double>(pars(348), pars(349));
   ZEL(0,1) = std::complex<double>(pars(350), pars(351));
   ZEL(0,2) = std::complex<double>(pars(352), pars(353));
   ZEL(1,0) = std::complex<double>(pars(354), pars(355));
   ZEL(1,1) = std::complex<double>(pars(356), pars(357));
   ZEL(1,2) = std::complex<double>(pars(358), pars(359));
   ZEL(2,0) = std::complex<double>(pars(360), pars(361));
   ZEL(2,1) = std::complex<double>(pars(362), pars(363));
   ZEL(2,2) = std::complex<double>(pars(364), pars(365));
   ZER(0,0) = std::complex<double>(pars(366), pars(367));
   ZER(0,1) = std::complex<double>(pars(368), pars(369));
   ZER(0,2) = std::complex<double>(pars(370), pars(371));
   ZER(1,0) = std::complex<double>(pars(372), pars(373));
   ZER(1,1) = std::complex<double>(pars(374), pars(375));
   ZER(1,2) = std::complex<double>(pars(376), pars(377));
   ZER(2,0) = std::complex<double>(pars(378), pars(379));
   ZER(2,1) = std::complex<double>(pars(380), pars(381));
   ZER(2,2) = std::complex<double>(pars(382), pars(383));
   ZDL(0,0) = std::complex<double>(pars(384), pars(385));
   ZDL(0,1) = std::complex<double>(pars(386), pars(387));
   ZDL(0,2) = std::complex<double>(pars(388), pars(389));
   ZDL(1,0) = std::complex<double>(pars(390), pars(391));
   ZDL(1,1) = std::complex<double>(pars(392), pars(393));
   ZDL(1,2) = std::complex<double>(pars(394), pars(395));
   ZDL(2,0) = std::complex<double>(pars(396), pars(397));
   ZDL(2,1) = std::complex<double>(pars(398), pars(399));
   ZDL(2,2) = std::complex<double>(pars(400), pars(401));
   ZDR(0,0) = std::complex<double>(pars(402), pars(403));
   ZDR(0,1) = std::complex<double>(pars(404), pars(405));
   ZDR(0,2) = std::complex<double>(pars(406), pars(407));
   ZDR(1,0) = std::complex<double>(pars(408), pars(409));
   ZDR(1,1) = std::complex<double>(pars(410), pars(411));
   ZDR(1,2) = std::complex<double>(pars(412), pars(413));
   ZDR(2,0) = std::complex<double>(pars(414), pars(415));
   ZDR(2,1) = std::complex<double>(pars(416), pars(417));
   ZDR(2,2) = std::complex<double>(pars(418), pars(419));
   ZUL(0,0) = std::complex<double>(pars(420), pars(421));
   ZUL(0,1) = std::complex<double>(pars(422), pars(423));
   ZUL(0,2) = std::complex<double>(pars(424), pars(425));
   ZUL(1,0) = std::complex<double>(pars(426), pars(427));
   ZUL(1,1) = std::complex<double>(pars(428), pars(429));
   ZUL(1,2) = std::complex<double>(pars(430), pars(431));
   ZUL(2,0) = std::complex<double>(pars(432), pars(433));
   ZUL(2,1) = std::complex<double>(pars(434), pars(435));
   ZUL(2,2) = std::complex<double>(pars(436), pars(437));
   ZUR(0,0) = std::complex<double>(pars(438), pars(439));
   ZUR(0,1) = std::complex<double>(pars(440), pars(441));
   ZUR(0,2) = std::complex<double>(pars(442), pars(443));
   ZUR(1,0) = std::complex<double>(pars(444), pars(445));
   ZUR(1,1) = std::complex<double>(pars(446), pars(447));
   ZUR(1,2) = std::complex<double>(pars(448), pars(449));
   ZUR(2,0) = std::complex<double>(pars(450), pars(451));
   ZUR(2,1) = std::complex<double>(pars(452), pars(453));
   ZUR(2,2) = std::complex<double>(pars(454), pars(455));
   ZZ(0,0) = pars(456);
   ZZ(0,1) = pars(457);
   ZZ(0,2) = pars(458);
   ZZ(1,0) = pars(459);
   ZZ(1,1) = pars(460);
   ZZ(1,2) = pars(461);
   ZZ(2,0) = pars(462);
   ZZ(2,1) = pars(463);
   ZZ(2,2) = pars(464);

}

Eigen::ArrayXd UMSSM_physical::get_masses() const
{
   Eigen::ArrayXd pars(58);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MSd(0);
   pars(3) = MSd(1);
   pars(4) = MSd(2);
   pars(5) = MSd(3);
   pars(6) = MSd(4);
   pars(7) = MSd(5);
   pars(8) = MSv(0);
   pars(9) = MSv(1);
   pars(10) = MSv(2);
   pars(11) = MSv(3);
   pars(12) = MSv(4);
   pars(13) = MSv(5);
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
   pars(34) = MChi(0);
   pars(35) = MChi(1);
   pars(36) = MChi(2);
   pars(37) = MChi(3);
   pars(38) = MChi(4);
   pars(39) = MChi(5);
   pars(40) = MFv(0);
   pars(41) = MFv(1);
   pars(42) = MFv(2);
   pars(43) = MCha(0);
   pars(44) = MCha(1);
   pars(45) = MFe(0);
   pars(46) = MFe(1);
   pars(47) = MFe(2);
   pars(48) = MFd(0);
   pars(49) = MFd(1);
   pars(50) = MFd(2);
   pars(51) = MFu(0);
   pars(52) = MFu(1);
   pars(53) = MFu(2);
   pars(54) = MVWm;
   pars(55) = MVP;
   pars(56) = MVZ;
   pars(57) = MVZp;

   return pars;
}

void UMSSM_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MSd(0) = pars(2);
   MSd(1) = pars(3);
   MSd(2) = pars(4);
   MSd(3) = pars(5);
   MSd(4) = pars(6);
   MSd(5) = pars(7);
   MSv(0) = pars(8);
   MSv(1) = pars(9);
   MSv(2) = pars(10);
   MSv(3) = pars(11);
   MSv(4) = pars(12);
   MSv(5) = pars(13);
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
   MChi(0) = pars(34);
   MChi(1) = pars(35);
   MChi(2) = pars(36);
   MChi(3) = pars(37);
   MChi(4) = pars(38);
   MChi(5) = pars(39);
   MFv(0) = pars(40);
   MFv(1) = pars(41);
   MFv(2) = pars(42);
   MCha(0) = pars(43);
   MCha(1) = pars(44);
   MFe(0) = pars(45);
   MFe(1) = pars(46);
   MFe(2) = pars(47);
   MFd(0) = pars(48);
   MFd(1) = pars(49);
   MFd(2) = pars(50);
   MFu(0) = pars(51);
   MFu(1) = pars(52);
   MFu(2) = pars(53);
   MVWm = pars(54);
   MVP = pars(55);
   MVZ = pars(56);
   MVZp = pars(57);

}

void UMSSM_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
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
   ostr << "ZVL = " << ZVL << '\n';
   ostr << "ZVR = " << ZVR << '\n';
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

std::ostream& operator<<(std::ostream& ostr, const UMSSM_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
