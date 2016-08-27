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

// File generated at Sat 27 Aug 2016 12:41:05

#include "MRSSM_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

MRSSM_physical::MRSSM_physical()
   :
    MVG(0), MGlu(0), MFv(Eigen::Array<double,3,1>::Zero()), MSRdp(0), MSRum(0)
       , MsigmaO(0), MphiO(0), MSd(Eigen::Array<double,6,1>::Zero()), MSv(
       Eigen::Array<double,3,1>::Zero()), MSu(Eigen::Array<double,6,1>::Zero()),
       MSe(Eigen::Array<double,6,1>::Zero()), Mhh(Eigen::Array<double,4,1>::Zero()
       ), MAh(Eigen::Array<double,4,1>::Zero()), MRh(Eigen::Array<double,2,1>
       ::Zero()), MHpm(Eigen::Array<double,4,1>::Zero()), MChi(Eigen::Array<double
       ,4,1>::Zero()), MCha1(Eigen::Array<double,2,1>::Zero()), MCha2(Eigen::Array
       <double,2,1>::Zero()), MFe(Eigen::Array<double,3,1>::Zero()), MFd(
       Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<double,3,1>::Zero()),
       MVWm(0), MVP(0), MVZ(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,3,3>::Zero(
      )), ZU(Eigen::Matrix<double,6,6>::Zero()), ZE(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,4,4>::Zero()), ZA(Eigen::Matrix<double,4,
      4>::Zero()), ZHR(Eigen::Matrix<double,2,2>::Zero()), ZP(Eigen::Matrix<double
      ,4,4>::Zero()), ZN1(Eigen::Matrix<std::complex<double>,4,4>::Zero()), ZN2(
      Eigen::Matrix<std::complex<double>,4,4>::Zero()), UM1(Eigen::Matrix<
      std::complex<double>,2,2>::Zero()), UP1(Eigen::Matrix<std::complex<double>,2
      ,2>::Zero()), UM2(Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP2(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), ZEL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZER(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZDL(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDR(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZUR(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())

{
}

void MRSSM_physical::clear()
{
   MVG = 0.;
   MGlu = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   MSRdp = 0.;
   MSRum = 0.;
   MsigmaO = 0.;
   MphiO = 0.;
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,3,1>::Zero();
   ZV = Eigen::Matrix<double,3,3>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,4,1>::Zero();
   ZH = Eigen::Matrix<double,4,4>::Zero();
   MAh = Eigen::Matrix<double,4,1>::Zero();
   ZA = Eigen::Matrix<double,4,4>::Zero();
   MRh = Eigen::Matrix<double,2,1>::Zero();
   ZHR = Eigen::Matrix<double,2,2>::Zero();
   MHpm = Eigen::Matrix<double,4,1>::Zero();
   ZP = Eigen::Matrix<double,4,4>::Zero();
   MChi = Eigen::Matrix<double,4,1>::Zero();
   ZN1 = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   ZN2 = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MCha1 = Eigen::Matrix<double,2,1>::Zero();
   UM1 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP1 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   MCha2 = Eigen::Matrix<double,2,1>::Zero();
   UM2 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
   UP2 = Eigen::Matrix<std::complex<double>,2,2>::Zero();
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
void MRSSM_physical::convert_to_hk()
{

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void MRSSM_physical::convert_to_slha()
{

}

Eigen::ArrayXd MRSSM_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(441);

   pars(64) = ZD(0,0);
   pars(65) = ZD(0,1);
   pars(66) = ZD(0,2);
   pars(67) = ZD(0,3);
   pars(68) = ZD(0,4);
   pars(69) = ZD(0,5);
   pars(70) = ZD(1,0);
   pars(71) = ZD(1,1);
   pars(72) = ZD(1,2);
   pars(73) = ZD(1,3);
   pars(74) = ZD(1,4);
   pars(75) = ZD(1,5);
   pars(76) = ZD(2,0);
   pars(77) = ZD(2,1);
   pars(78) = ZD(2,2);
   pars(79) = ZD(2,3);
   pars(80) = ZD(2,4);
   pars(81) = ZD(2,5);
   pars(82) = ZD(3,0);
   pars(83) = ZD(3,1);
   pars(84) = ZD(3,2);
   pars(85) = ZD(3,3);
   pars(86) = ZD(3,4);
   pars(87) = ZD(3,5);
   pars(88) = ZD(4,0);
   pars(89) = ZD(4,1);
   pars(90) = ZD(4,2);
   pars(91) = ZD(4,3);
   pars(92) = ZD(4,4);
   pars(93) = ZD(4,5);
   pars(94) = ZD(5,0);
   pars(95) = ZD(5,1);
   pars(96) = ZD(5,2);
   pars(97) = ZD(5,3);
   pars(98) = ZD(5,4);
   pars(99) = ZD(5,5);
   pars(100) = ZV(0,0);
   pars(101) = ZV(0,1);
   pars(102) = ZV(0,2);
   pars(103) = ZV(1,0);
   pars(104) = ZV(1,1);
   pars(105) = ZV(1,2);
   pars(106) = ZV(2,0);
   pars(107) = ZV(2,1);
   pars(108) = ZV(2,2);
   pars(109) = ZU(0,0);
   pars(110) = ZU(0,1);
   pars(111) = ZU(0,2);
   pars(112) = ZU(0,3);
   pars(113) = ZU(0,4);
   pars(114) = ZU(0,5);
   pars(115) = ZU(1,0);
   pars(116) = ZU(1,1);
   pars(117) = ZU(1,2);
   pars(118) = ZU(1,3);
   pars(119) = ZU(1,4);
   pars(120) = ZU(1,5);
   pars(121) = ZU(2,0);
   pars(122) = ZU(2,1);
   pars(123) = ZU(2,2);
   pars(124) = ZU(2,3);
   pars(125) = ZU(2,4);
   pars(126) = ZU(2,5);
   pars(127) = ZU(3,0);
   pars(128) = ZU(3,1);
   pars(129) = ZU(3,2);
   pars(130) = ZU(3,3);
   pars(131) = ZU(3,4);
   pars(132) = ZU(3,5);
   pars(133) = ZU(4,0);
   pars(134) = ZU(4,1);
   pars(135) = ZU(4,2);
   pars(136) = ZU(4,3);
   pars(137) = ZU(4,4);
   pars(138) = ZU(4,5);
   pars(139) = ZU(5,0);
   pars(140) = ZU(5,1);
   pars(141) = ZU(5,2);
   pars(142) = ZU(5,3);
   pars(143) = ZU(5,4);
   pars(144) = ZU(5,5);
   pars(145) = ZE(0,0);
   pars(146) = ZE(0,1);
   pars(147) = ZE(0,2);
   pars(148) = ZE(0,3);
   pars(149) = ZE(0,4);
   pars(150) = ZE(0,5);
   pars(151) = ZE(1,0);
   pars(152) = ZE(1,1);
   pars(153) = ZE(1,2);
   pars(154) = ZE(1,3);
   pars(155) = ZE(1,4);
   pars(156) = ZE(1,5);
   pars(157) = ZE(2,0);
   pars(158) = ZE(2,1);
   pars(159) = ZE(2,2);
   pars(160) = ZE(2,3);
   pars(161) = ZE(2,4);
   pars(162) = ZE(2,5);
   pars(163) = ZE(3,0);
   pars(164) = ZE(3,1);
   pars(165) = ZE(3,2);
   pars(166) = ZE(3,3);
   pars(167) = ZE(3,4);
   pars(168) = ZE(3,5);
   pars(169) = ZE(4,0);
   pars(170) = ZE(4,1);
   pars(171) = ZE(4,2);
   pars(172) = ZE(4,3);
   pars(173) = ZE(4,4);
   pars(174) = ZE(4,5);
   pars(175) = ZE(5,0);
   pars(176) = ZE(5,1);
   pars(177) = ZE(5,2);
   pars(178) = ZE(5,3);
   pars(179) = ZE(5,4);
   pars(180) = ZE(5,5);
   pars(181) = ZH(0,0);
   pars(182) = ZH(0,1);
   pars(183) = ZH(0,2);
   pars(184) = ZH(0,3);
   pars(185) = ZH(1,0);
   pars(186) = ZH(1,1);
   pars(187) = ZH(1,2);
   pars(188) = ZH(1,3);
   pars(189) = ZH(2,0);
   pars(190) = ZH(2,1);
   pars(191) = ZH(2,2);
   pars(192) = ZH(2,3);
   pars(193) = ZH(3,0);
   pars(194) = ZH(3,1);
   pars(195) = ZH(3,2);
   pars(196) = ZH(3,3);
   pars(197) = ZA(0,0);
   pars(198) = ZA(0,1);
   pars(199) = ZA(0,2);
   pars(200) = ZA(0,3);
   pars(201) = ZA(1,0);
   pars(202) = ZA(1,1);
   pars(203) = ZA(1,2);
   pars(204) = ZA(1,3);
   pars(205) = ZA(2,0);
   pars(206) = ZA(2,1);
   pars(207) = ZA(2,2);
   pars(208) = ZA(2,3);
   pars(209) = ZA(3,0);
   pars(210) = ZA(3,1);
   pars(211) = ZA(3,2);
   pars(212) = ZA(3,3);
   pars(213) = ZHR(0,0);
   pars(214) = ZHR(0,1);
   pars(215) = ZHR(1,0);
   pars(216) = ZHR(1,1);
   pars(217) = ZP(0,0);
   pars(218) = ZP(0,1);
   pars(219) = ZP(0,2);
   pars(220) = ZP(0,3);
   pars(221) = ZP(1,0);
   pars(222) = ZP(1,1);
   pars(223) = ZP(1,2);
   pars(224) = ZP(1,3);
   pars(225) = ZP(2,0);
   pars(226) = ZP(2,1);
   pars(227) = ZP(2,2);
   pars(228) = ZP(2,3);
   pars(229) = ZP(3,0);
   pars(230) = ZP(3,1);
   pars(231) = ZP(3,2);
   pars(232) = ZP(3,3);
   pars(233) = Re(ZN1(0,0));
   pars(234) = Im(ZN1(0,0));
   pars(235) = Re(ZN1(0,1));
   pars(236) = Im(ZN1(0,1));
   pars(237) = Re(ZN1(0,2));
   pars(238) = Im(ZN1(0,2));
   pars(239) = Re(ZN1(0,3));
   pars(240) = Im(ZN1(0,3));
   pars(241) = Re(ZN1(1,0));
   pars(242) = Im(ZN1(1,0));
   pars(243) = Re(ZN1(1,1));
   pars(244) = Im(ZN1(1,1));
   pars(245) = Re(ZN1(1,2));
   pars(246) = Im(ZN1(1,2));
   pars(247) = Re(ZN1(1,3));
   pars(248) = Im(ZN1(1,3));
   pars(249) = Re(ZN1(2,0));
   pars(250) = Im(ZN1(2,0));
   pars(251) = Re(ZN1(2,1));
   pars(252) = Im(ZN1(2,1));
   pars(253) = Re(ZN1(2,2));
   pars(254) = Im(ZN1(2,2));
   pars(255) = Re(ZN1(2,3));
   pars(256) = Im(ZN1(2,3));
   pars(257) = Re(ZN1(3,0));
   pars(258) = Im(ZN1(3,0));
   pars(259) = Re(ZN1(3,1));
   pars(260) = Im(ZN1(3,1));
   pars(261) = Re(ZN1(3,2));
   pars(262) = Im(ZN1(3,2));
   pars(263) = Re(ZN1(3,3));
   pars(264) = Im(ZN1(3,3));
   pars(265) = Re(ZN2(0,0));
   pars(266) = Im(ZN2(0,0));
   pars(267) = Re(ZN2(0,1));
   pars(268) = Im(ZN2(0,1));
   pars(269) = Re(ZN2(0,2));
   pars(270) = Im(ZN2(0,2));
   pars(271) = Re(ZN2(0,3));
   pars(272) = Im(ZN2(0,3));
   pars(273) = Re(ZN2(1,0));
   pars(274) = Im(ZN2(1,0));
   pars(275) = Re(ZN2(1,1));
   pars(276) = Im(ZN2(1,1));
   pars(277) = Re(ZN2(1,2));
   pars(278) = Im(ZN2(1,2));
   pars(279) = Re(ZN2(1,3));
   pars(280) = Im(ZN2(1,3));
   pars(281) = Re(ZN2(2,0));
   pars(282) = Im(ZN2(2,0));
   pars(283) = Re(ZN2(2,1));
   pars(284) = Im(ZN2(2,1));
   pars(285) = Re(ZN2(2,2));
   pars(286) = Im(ZN2(2,2));
   pars(287) = Re(ZN2(2,3));
   pars(288) = Im(ZN2(2,3));
   pars(289) = Re(ZN2(3,0));
   pars(290) = Im(ZN2(3,0));
   pars(291) = Re(ZN2(3,1));
   pars(292) = Im(ZN2(3,1));
   pars(293) = Re(ZN2(3,2));
   pars(294) = Im(ZN2(3,2));
   pars(295) = Re(ZN2(3,3));
   pars(296) = Im(ZN2(3,3));
   pars(297) = Re(UM1(0,0));
   pars(298) = Im(UM1(0,0));
   pars(299) = Re(UM1(0,1));
   pars(300) = Im(UM1(0,1));
   pars(301) = Re(UM1(1,0));
   pars(302) = Im(UM1(1,0));
   pars(303) = Re(UM1(1,1));
   pars(304) = Im(UM1(1,1));
   pars(305) = Re(UP1(0,0));
   pars(306) = Im(UP1(0,0));
   pars(307) = Re(UP1(0,1));
   pars(308) = Im(UP1(0,1));
   pars(309) = Re(UP1(1,0));
   pars(310) = Im(UP1(1,0));
   pars(311) = Re(UP1(1,1));
   pars(312) = Im(UP1(1,1));
   pars(313) = Re(UM2(0,0));
   pars(314) = Im(UM2(0,0));
   pars(315) = Re(UM2(0,1));
   pars(316) = Im(UM2(0,1));
   pars(317) = Re(UM2(1,0));
   pars(318) = Im(UM2(1,0));
   pars(319) = Re(UM2(1,1));
   pars(320) = Im(UM2(1,1));
   pars(321) = Re(UP2(0,0));
   pars(322) = Im(UP2(0,0));
   pars(323) = Re(UP2(0,1));
   pars(324) = Im(UP2(0,1));
   pars(325) = Re(UP2(1,0));
   pars(326) = Im(UP2(1,0));
   pars(327) = Re(UP2(1,1));
   pars(328) = Im(UP2(1,1));
   pars(329) = Re(ZEL(0,0));
   pars(330) = Im(ZEL(0,0));
   pars(331) = Re(ZEL(0,1));
   pars(332) = Im(ZEL(0,1));
   pars(333) = Re(ZEL(0,2));
   pars(334) = Im(ZEL(0,2));
   pars(335) = Re(ZEL(1,0));
   pars(336) = Im(ZEL(1,0));
   pars(337) = Re(ZEL(1,1));
   pars(338) = Im(ZEL(1,1));
   pars(339) = Re(ZEL(1,2));
   pars(340) = Im(ZEL(1,2));
   pars(341) = Re(ZEL(2,0));
   pars(342) = Im(ZEL(2,0));
   pars(343) = Re(ZEL(2,1));
   pars(344) = Im(ZEL(2,1));
   pars(345) = Re(ZEL(2,2));
   pars(346) = Im(ZEL(2,2));
   pars(347) = Re(ZER(0,0));
   pars(348) = Im(ZER(0,0));
   pars(349) = Re(ZER(0,1));
   pars(350) = Im(ZER(0,1));
   pars(351) = Re(ZER(0,2));
   pars(352) = Im(ZER(0,2));
   pars(353) = Re(ZER(1,0));
   pars(354) = Im(ZER(1,0));
   pars(355) = Re(ZER(1,1));
   pars(356) = Im(ZER(1,1));
   pars(357) = Re(ZER(1,2));
   pars(358) = Im(ZER(1,2));
   pars(359) = Re(ZER(2,0));
   pars(360) = Im(ZER(2,0));
   pars(361) = Re(ZER(2,1));
   pars(362) = Im(ZER(2,1));
   pars(363) = Re(ZER(2,2));
   pars(364) = Im(ZER(2,2));
   pars(365) = Re(ZDL(0,0));
   pars(366) = Im(ZDL(0,0));
   pars(367) = Re(ZDL(0,1));
   pars(368) = Im(ZDL(0,1));
   pars(369) = Re(ZDL(0,2));
   pars(370) = Im(ZDL(0,2));
   pars(371) = Re(ZDL(1,0));
   pars(372) = Im(ZDL(1,0));
   pars(373) = Re(ZDL(1,1));
   pars(374) = Im(ZDL(1,1));
   pars(375) = Re(ZDL(1,2));
   pars(376) = Im(ZDL(1,2));
   pars(377) = Re(ZDL(2,0));
   pars(378) = Im(ZDL(2,0));
   pars(379) = Re(ZDL(2,1));
   pars(380) = Im(ZDL(2,1));
   pars(381) = Re(ZDL(2,2));
   pars(382) = Im(ZDL(2,2));
   pars(383) = Re(ZDR(0,0));
   pars(384) = Im(ZDR(0,0));
   pars(385) = Re(ZDR(0,1));
   pars(386) = Im(ZDR(0,1));
   pars(387) = Re(ZDR(0,2));
   pars(388) = Im(ZDR(0,2));
   pars(389) = Re(ZDR(1,0));
   pars(390) = Im(ZDR(1,0));
   pars(391) = Re(ZDR(1,1));
   pars(392) = Im(ZDR(1,1));
   pars(393) = Re(ZDR(1,2));
   pars(394) = Im(ZDR(1,2));
   pars(395) = Re(ZDR(2,0));
   pars(396) = Im(ZDR(2,0));
   pars(397) = Re(ZDR(2,1));
   pars(398) = Im(ZDR(2,1));
   pars(399) = Re(ZDR(2,2));
   pars(400) = Im(ZDR(2,2));
   pars(401) = Re(ZUL(0,0));
   pars(402) = Im(ZUL(0,0));
   pars(403) = Re(ZUL(0,1));
   pars(404) = Im(ZUL(0,1));
   pars(405) = Re(ZUL(0,2));
   pars(406) = Im(ZUL(0,2));
   pars(407) = Re(ZUL(1,0));
   pars(408) = Im(ZUL(1,0));
   pars(409) = Re(ZUL(1,1));
   pars(410) = Im(ZUL(1,1));
   pars(411) = Re(ZUL(1,2));
   pars(412) = Im(ZUL(1,2));
   pars(413) = Re(ZUL(2,0));
   pars(414) = Im(ZUL(2,0));
   pars(415) = Re(ZUL(2,1));
   pars(416) = Im(ZUL(2,1));
   pars(417) = Re(ZUL(2,2));
   pars(418) = Im(ZUL(2,2));
   pars(419) = Re(ZUR(0,0));
   pars(420) = Im(ZUR(0,0));
   pars(421) = Re(ZUR(0,1));
   pars(422) = Im(ZUR(0,1));
   pars(423) = Re(ZUR(0,2));
   pars(424) = Im(ZUR(0,2));
   pars(425) = Re(ZUR(1,0));
   pars(426) = Im(ZUR(1,0));
   pars(427) = Re(ZUR(1,1));
   pars(428) = Im(ZUR(1,1));
   pars(429) = Re(ZUR(1,2));
   pars(430) = Im(ZUR(1,2));
   pars(431) = Re(ZUR(2,0));
   pars(432) = Im(ZUR(2,0));
   pars(433) = Re(ZUR(2,1));
   pars(434) = Im(ZUR(2,1));
   pars(435) = Re(ZUR(2,2));
   pars(436) = Im(ZUR(2,2));
   pars(437) = ZZ(0,0);
   pars(438) = ZZ(0,1);
   pars(439) = ZZ(1,0);
   pars(440) = ZZ(1,1);


   return pars;
}

void MRSSM_physical::set(const Eigen::ArrayXd& pars)
{
   set_masses(pars);

   ZD(0,0) = pars(64);
   ZD(0,1) = pars(65);
   ZD(0,2) = pars(66);
   ZD(0,3) = pars(67);
   ZD(0,4) = pars(68);
   ZD(0,5) = pars(69);
   ZD(1,0) = pars(70);
   ZD(1,1) = pars(71);
   ZD(1,2) = pars(72);
   ZD(1,3) = pars(73);
   ZD(1,4) = pars(74);
   ZD(1,5) = pars(75);
   ZD(2,0) = pars(76);
   ZD(2,1) = pars(77);
   ZD(2,2) = pars(78);
   ZD(2,3) = pars(79);
   ZD(2,4) = pars(80);
   ZD(2,5) = pars(81);
   ZD(3,0) = pars(82);
   ZD(3,1) = pars(83);
   ZD(3,2) = pars(84);
   ZD(3,3) = pars(85);
   ZD(3,4) = pars(86);
   ZD(3,5) = pars(87);
   ZD(4,0) = pars(88);
   ZD(4,1) = pars(89);
   ZD(4,2) = pars(90);
   ZD(4,3) = pars(91);
   ZD(4,4) = pars(92);
   ZD(4,5) = pars(93);
   ZD(5,0) = pars(94);
   ZD(5,1) = pars(95);
   ZD(5,2) = pars(96);
   ZD(5,3) = pars(97);
   ZD(5,4) = pars(98);
   ZD(5,5) = pars(99);
   ZV(0,0) = pars(100);
   ZV(0,1) = pars(101);
   ZV(0,2) = pars(102);
   ZV(1,0) = pars(103);
   ZV(1,1) = pars(104);
   ZV(1,2) = pars(105);
   ZV(2,0) = pars(106);
   ZV(2,1) = pars(107);
   ZV(2,2) = pars(108);
   ZU(0,0) = pars(109);
   ZU(0,1) = pars(110);
   ZU(0,2) = pars(111);
   ZU(0,3) = pars(112);
   ZU(0,4) = pars(113);
   ZU(0,5) = pars(114);
   ZU(1,0) = pars(115);
   ZU(1,1) = pars(116);
   ZU(1,2) = pars(117);
   ZU(1,3) = pars(118);
   ZU(1,4) = pars(119);
   ZU(1,5) = pars(120);
   ZU(2,0) = pars(121);
   ZU(2,1) = pars(122);
   ZU(2,2) = pars(123);
   ZU(2,3) = pars(124);
   ZU(2,4) = pars(125);
   ZU(2,5) = pars(126);
   ZU(3,0) = pars(127);
   ZU(3,1) = pars(128);
   ZU(3,2) = pars(129);
   ZU(3,3) = pars(130);
   ZU(3,4) = pars(131);
   ZU(3,5) = pars(132);
   ZU(4,0) = pars(133);
   ZU(4,1) = pars(134);
   ZU(4,2) = pars(135);
   ZU(4,3) = pars(136);
   ZU(4,4) = pars(137);
   ZU(4,5) = pars(138);
   ZU(5,0) = pars(139);
   ZU(5,1) = pars(140);
   ZU(5,2) = pars(141);
   ZU(5,3) = pars(142);
   ZU(5,4) = pars(143);
   ZU(5,5) = pars(144);
   ZE(0,0) = pars(145);
   ZE(0,1) = pars(146);
   ZE(0,2) = pars(147);
   ZE(0,3) = pars(148);
   ZE(0,4) = pars(149);
   ZE(0,5) = pars(150);
   ZE(1,0) = pars(151);
   ZE(1,1) = pars(152);
   ZE(1,2) = pars(153);
   ZE(1,3) = pars(154);
   ZE(1,4) = pars(155);
   ZE(1,5) = pars(156);
   ZE(2,0) = pars(157);
   ZE(2,1) = pars(158);
   ZE(2,2) = pars(159);
   ZE(2,3) = pars(160);
   ZE(2,4) = pars(161);
   ZE(2,5) = pars(162);
   ZE(3,0) = pars(163);
   ZE(3,1) = pars(164);
   ZE(3,2) = pars(165);
   ZE(3,3) = pars(166);
   ZE(3,4) = pars(167);
   ZE(3,5) = pars(168);
   ZE(4,0) = pars(169);
   ZE(4,1) = pars(170);
   ZE(4,2) = pars(171);
   ZE(4,3) = pars(172);
   ZE(4,4) = pars(173);
   ZE(4,5) = pars(174);
   ZE(5,0) = pars(175);
   ZE(5,1) = pars(176);
   ZE(5,2) = pars(177);
   ZE(5,3) = pars(178);
   ZE(5,4) = pars(179);
   ZE(5,5) = pars(180);
   ZH(0,0) = pars(181);
   ZH(0,1) = pars(182);
   ZH(0,2) = pars(183);
   ZH(0,3) = pars(184);
   ZH(1,0) = pars(185);
   ZH(1,1) = pars(186);
   ZH(1,2) = pars(187);
   ZH(1,3) = pars(188);
   ZH(2,0) = pars(189);
   ZH(2,1) = pars(190);
   ZH(2,2) = pars(191);
   ZH(2,3) = pars(192);
   ZH(3,0) = pars(193);
   ZH(3,1) = pars(194);
   ZH(3,2) = pars(195);
   ZH(3,3) = pars(196);
   ZA(0,0) = pars(197);
   ZA(0,1) = pars(198);
   ZA(0,2) = pars(199);
   ZA(0,3) = pars(200);
   ZA(1,0) = pars(201);
   ZA(1,1) = pars(202);
   ZA(1,2) = pars(203);
   ZA(1,3) = pars(204);
   ZA(2,0) = pars(205);
   ZA(2,1) = pars(206);
   ZA(2,2) = pars(207);
   ZA(2,3) = pars(208);
   ZA(3,0) = pars(209);
   ZA(3,1) = pars(210);
   ZA(3,2) = pars(211);
   ZA(3,3) = pars(212);
   ZHR(0,0) = pars(213);
   ZHR(0,1) = pars(214);
   ZHR(1,0) = pars(215);
   ZHR(1,1) = pars(216);
   ZP(0,0) = pars(217);
   ZP(0,1) = pars(218);
   ZP(0,2) = pars(219);
   ZP(0,3) = pars(220);
   ZP(1,0) = pars(221);
   ZP(1,1) = pars(222);
   ZP(1,2) = pars(223);
   ZP(1,3) = pars(224);
   ZP(2,0) = pars(225);
   ZP(2,1) = pars(226);
   ZP(2,2) = pars(227);
   ZP(2,3) = pars(228);
   ZP(3,0) = pars(229);
   ZP(3,1) = pars(230);
   ZP(3,2) = pars(231);
   ZP(3,3) = pars(232);
   ZN1(0,0) = std::complex<double>(pars(233), pars(234));
   ZN1(0,1) = std::complex<double>(pars(235), pars(236));
   ZN1(0,2) = std::complex<double>(pars(237), pars(238));
   ZN1(0,3) = std::complex<double>(pars(239), pars(240));
   ZN1(1,0) = std::complex<double>(pars(241), pars(242));
   ZN1(1,1) = std::complex<double>(pars(243), pars(244));
   ZN1(1,2) = std::complex<double>(pars(245), pars(246));
   ZN1(1,3) = std::complex<double>(pars(247), pars(248));
   ZN1(2,0) = std::complex<double>(pars(249), pars(250));
   ZN1(2,1) = std::complex<double>(pars(251), pars(252));
   ZN1(2,2) = std::complex<double>(pars(253), pars(254));
   ZN1(2,3) = std::complex<double>(pars(255), pars(256));
   ZN1(3,0) = std::complex<double>(pars(257), pars(258));
   ZN1(3,1) = std::complex<double>(pars(259), pars(260));
   ZN1(3,2) = std::complex<double>(pars(261), pars(262));
   ZN1(3,3) = std::complex<double>(pars(263), pars(264));
   ZN2(0,0) = std::complex<double>(pars(265), pars(266));
   ZN2(0,1) = std::complex<double>(pars(267), pars(268));
   ZN2(0,2) = std::complex<double>(pars(269), pars(270));
   ZN2(0,3) = std::complex<double>(pars(271), pars(272));
   ZN2(1,0) = std::complex<double>(pars(273), pars(274));
   ZN2(1,1) = std::complex<double>(pars(275), pars(276));
   ZN2(1,2) = std::complex<double>(pars(277), pars(278));
   ZN2(1,3) = std::complex<double>(pars(279), pars(280));
   ZN2(2,0) = std::complex<double>(pars(281), pars(282));
   ZN2(2,1) = std::complex<double>(pars(283), pars(284));
   ZN2(2,2) = std::complex<double>(pars(285), pars(286));
   ZN2(2,3) = std::complex<double>(pars(287), pars(288));
   ZN2(3,0) = std::complex<double>(pars(289), pars(290));
   ZN2(3,1) = std::complex<double>(pars(291), pars(292));
   ZN2(3,2) = std::complex<double>(pars(293), pars(294));
   ZN2(3,3) = std::complex<double>(pars(295), pars(296));
   UM1(0,0) = std::complex<double>(pars(297), pars(298));
   UM1(0,1) = std::complex<double>(pars(299), pars(300));
   UM1(1,0) = std::complex<double>(pars(301), pars(302));
   UM1(1,1) = std::complex<double>(pars(303), pars(304));
   UP1(0,0) = std::complex<double>(pars(305), pars(306));
   UP1(0,1) = std::complex<double>(pars(307), pars(308));
   UP1(1,0) = std::complex<double>(pars(309), pars(310));
   UP1(1,1) = std::complex<double>(pars(311), pars(312));
   UM2(0,0) = std::complex<double>(pars(313), pars(314));
   UM2(0,1) = std::complex<double>(pars(315), pars(316));
   UM2(1,0) = std::complex<double>(pars(317), pars(318));
   UM2(1,1) = std::complex<double>(pars(319), pars(320));
   UP2(0,0) = std::complex<double>(pars(321), pars(322));
   UP2(0,1) = std::complex<double>(pars(323), pars(324));
   UP2(1,0) = std::complex<double>(pars(325), pars(326));
   UP2(1,1) = std::complex<double>(pars(327), pars(328));
   ZEL(0,0) = std::complex<double>(pars(329), pars(330));
   ZEL(0,1) = std::complex<double>(pars(331), pars(332));
   ZEL(0,2) = std::complex<double>(pars(333), pars(334));
   ZEL(1,0) = std::complex<double>(pars(335), pars(336));
   ZEL(1,1) = std::complex<double>(pars(337), pars(338));
   ZEL(1,2) = std::complex<double>(pars(339), pars(340));
   ZEL(2,0) = std::complex<double>(pars(341), pars(342));
   ZEL(2,1) = std::complex<double>(pars(343), pars(344));
   ZEL(2,2) = std::complex<double>(pars(345), pars(346));
   ZER(0,0) = std::complex<double>(pars(347), pars(348));
   ZER(0,1) = std::complex<double>(pars(349), pars(350));
   ZER(0,2) = std::complex<double>(pars(351), pars(352));
   ZER(1,0) = std::complex<double>(pars(353), pars(354));
   ZER(1,1) = std::complex<double>(pars(355), pars(356));
   ZER(1,2) = std::complex<double>(pars(357), pars(358));
   ZER(2,0) = std::complex<double>(pars(359), pars(360));
   ZER(2,1) = std::complex<double>(pars(361), pars(362));
   ZER(2,2) = std::complex<double>(pars(363), pars(364));
   ZDL(0,0) = std::complex<double>(pars(365), pars(366));
   ZDL(0,1) = std::complex<double>(pars(367), pars(368));
   ZDL(0,2) = std::complex<double>(pars(369), pars(370));
   ZDL(1,0) = std::complex<double>(pars(371), pars(372));
   ZDL(1,1) = std::complex<double>(pars(373), pars(374));
   ZDL(1,2) = std::complex<double>(pars(375), pars(376));
   ZDL(2,0) = std::complex<double>(pars(377), pars(378));
   ZDL(2,1) = std::complex<double>(pars(379), pars(380));
   ZDL(2,2) = std::complex<double>(pars(381), pars(382));
   ZDR(0,0) = std::complex<double>(pars(383), pars(384));
   ZDR(0,1) = std::complex<double>(pars(385), pars(386));
   ZDR(0,2) = std::complex<double>(pars(387), pars(388));
   ZDR(1,0) = std::complex<double>(pars(389), pars(390));
   ZDR(1,1) = std::complex<double>(pars(391), pars(392));
   ZDR(1,2) = std::complex<double>(pars(393), pars(394));
   ZDR(2,0) = std::complex<double>(pars(395), pars(396));
   ZDR(2,1) = std::complex<double>(pars(397), pars(398));
   ZDR(2,2) = std::complex<double>(pars(399), pars(400));
   ZUL(0,0) = std::complex<double>(pars(401), pars(402));
   ZUL(0,1) = std::complex<double>(pars(403), pars(404));
   ZUL(0,2) = std::complex<double>(pars(405), pars(406));
   ZUL(1,0) = std::complex<double>(pars(407), pars(408));
   ZUL(1,1) = std::complex<double>(pars(409), pars(410));
   ZUL(1,2) = std::complex<double>(pars(411), pars(412));
   ZUL(2,0) = std::complex<double>(pars(413), pars(414));
   ZUL(2,1) = std::complex<double>(pars(415), pars(416));
   ZUL(2,2) = std::complex<double>(pars(417), pars(418));
   ZUR(0,0) = std::complex<double>(pars(419), pars(420));
   ZUR(0,1) = std::complex<double>(pars(421), pars(422));
   ZUR(0,2) = std::complex<double>(pars(423), pars(424));
   ZUR(1,0) = std::complex<double>(pars(425), pars(426));
   ZUR(1,1) = std::complex<double>(pars(427), pars(428));
   ZUR(1,2) = std::complex<double>(pars(429), pars(430));
   ZUR(2,0) = std::complex<double>(pars(431), pars(432));
   ZUR(2,1) = std::complex<double>(pars(433), pars(434));
   ZUR(2,2) = std::complex<double>(pars(435), pars(436));
   ZZ(0,0) = pars(437);
   ZZ(0,1) = pars(438);
   ZZ(1,0) = pars(439);
   ZZ(1,1) = pars(440);

}

Eigen::ArrayXd MRSSM_physical::get_masses() const
{
   Eigen::ArrayXd pars(64);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MFv(0);
   pars(3) = MFv(1);
   pars(4) = MFv(2);
   pars(5) = MSRdp;
   pars(6) = MSRum;
   pars(7) = MsigmaO;
   pars(8) = MphiO;
   pars(9) = MSd(0);
   pars(10) = MSd(1);
   pars(11) = MSd(2);
   pars(12) = MSd(3);
   pars(13) = MSd(4);
   pars(14) = MSd(5);
   pars(15) = MSv(0);
   pars(16) = MSv(1);
   pars(17) = MSv(2);
   pars(18) = MSu(0);
   pars(19) = MSu(1);
   pars(20) = MSu(2);
   pars(21) = MSu(3);
   pars(22) = MSu(4);
   pars(23) = MSu(5);
   pars(24) = MSe(0);
   pars(25) = MSe(1);
   pars(26) = MSe(2);
   pars(27) = MSe(3);
   pars(28) = MSe(4);
   pars(29) = MSe(5);
   pars(30) = Mhh(0);
   pars(31) = Mhh(1);
   pars(32) = Mhh(2);
   pars(33) = Mhh(3);
   pars(34) = MAh(0);
   pars(35) = MAh(1);
   pars(36) = MAh(2);
   pars(37) = MAh(3);
   pars(38) = MRh(0);
   pars(39) = MRh(1);
   pars(40) = MHpm(0);
   pars(41) = MHpm(1);
   pars(42) = MHpm(2);
   pars(43) = MHpm(3);
   pars(44) = MChi(0);
   pars(45) = MChi(1);
   pars(46) = MChi(2);
   pars(47) = MChi(3);
   pars(48) = MCha1(0);
   pars(49) = MCha1(1);
   pars(50) = MCha2(0);
   pars(51) = MCha2(1);
   pars(52) = MFe(0);
   pars(53) = MFe(1);
   pars(54) = MFe(2);
   pars(55) = MFd(0);
   pars(56) = MFd(1);
   pars(57) = MFd(2);
   pars(58) = MFu(0);
   pars(59) = MFu(1);
   pars(60) = MFu(2);
   pars(61) = MVWm;
   pars(62) = MVP;
   pars(63) = MVZ;

   return pars;
}

void MRSSM_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MFv(0) = pars(2);
   MFv(1) = pars(3);
   MFv(2) = pars(4);
   MSRdp = pars(5);
   MSRum = pars(6);
   MsigmaO = pars(7);
   MphiO = pars(8);
   MSd(0) = pars(9);
   MSd(1) = pars(10);
   MSd(2) = pars(11);
   MSd(3) = pars(12);
   MSd(4) = pars(13);
   MSd(5) = pars(14);
   MSv(0) = pars(15);
   MSv(1) = pars(16);
   MSv(2) = pars(17);
   MSu(0) = pars(18);
   MSu(1) = pars(19);
   MSu(2) = pars(20);
   MSu(3) = pars(21);
   MSu(4) = pars(22);
   MSu(5) = pars(23);
   MSe(0) = pars(24);
   MSe(1) = pars(25);
   MSe(2) = pars(26);
   MSe(3) = pars(27);
   MSe(4) = pars(28);
   MSe(5) = pars(29);
   Mhh(0) = pars(30);
   Mhh(1) = pars(31);
   Mhh(2) = pars(32);
   Mhh(3) = pars(33);
   MAh(0) = pars(34);
   MAh(1) = pars(35);
   MAh(2) = pars(36);
   MAh(3) = pars(37);
   MRh(0) = pars(38);
   MRh(1) = pars(39);
   MHpm(0) = pars(40);
   MHpm(1) = pars(41);
   MHpm(2) = pars(42);
   MHpm(3) = pars(43);
   MChi(0) = pars(44);
   MChi(1) = pars(45);
   MChi(2) = pars(46);
   MChi(3) = pars(47);
   MCha1(0) = pars(48);
   MCha1(1) = pars(49);
   MCha2(0) = pars(50);
   MCha2(1) = pars(51);
   MFe(0) = pars(52);
   MFe(1) = pars(53);
   MFe(2) = pars(54);
   MFd(0) = pars(55);
   MFd(1) = pars(56);
   MFd(2) = pars(57);
   MFu(0) = pars(58);
   MFu(1) = pars(59);
   MFu(2) = pars(60);
   MVWm = pars(61);
   MVP = pars(62);
   MVZ = pars(63);

}

void MRSSM_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "MSRdp = " << MSRdp << '\n';
   ostr << "MSRum = " << MSRum << '\n';
   ostr << "MsigmaO = " << MsigmaO << '\n';
   ostr << "MphiO = " << MphiO << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "Mhh = " << Mhh.transpose() << '\n';
   ostr << "MAh = " << MAh.transpose() << '\n';
   ostr << "MRh = " << MRh.transpose() << '\n';
   ostr << "MHpm = " << MHpm.transpose() << '\n';
   ostr << "MChi = " << MChi.transpose() << '\n';
   ostr << "MCha1 = " << MCha1.transpose() << '\n';
   ostr << "MCha2 = " << MCha2.transpose() << '\n';
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
   ostr << "ZHR = " << ZHR << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN1 = " << ZN1 << '\n';
   ostr << "ZN2 = " << ZN2 << '\n';
   ostr << "UM1 = " << UM1 << '\n';
   ostr << "UP1 = " << UP1 << '\n';
   ostr << "UM2 = " << UM2 << '\n';
   ostr << "UP2 = " << UP2 << '\n';
   ostr << "ZEL = " << ZEL << '\n';
   ostr << "ZER = " << ZER << '\n';
   ostr << "ZDL = " << ZDL << '\n';
   ostr << "ZDR = " << ZDR << '\n';
   ostr << "ZUL = " << ZUL << '\n';
   ostr << "ZUR = " << ZUR << '\n';
   ostr << "ZZ = " << ZZ << '\n';

}

std::ostream& operator<<(std::ostream& ostr, const MRSSM_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
