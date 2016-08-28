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

// File generated at Sun 28 Aug 2016 15:18:29

#include "MSSMRHN_physical.hpp"
#include "slha_io.hpp"

#include <iostream>

#define LOCALPHYSICAL(p) p

namespace flexiblesusy {

MSSMRHN_physical::MSSMRHN_physical()
   :
    MVG(0), MGlu(0), MSd(Eigen::Array<double,6,1>::Zero()), MSu(Eigen::Array<
       double,6,1>::Zero()), MSe(Eigen::Array<double,6,1>::Zero()), MSv(
       Eigen::Array<double,6,1>::Zero()), Mhh(Eigen::Array<double,2,1>::Zero()),
       MAh(Eigen::Array<double,2,1>::Zero()), MHpm(Eigen::Array<double,2,1>::Zero(
       )), MChi(Eigen::Array<double,4,1>::Zero()), MFv(Eigen::Array<double,6,1>
       ::Zero()), MCha(Eigen::Array<double,2,1>::Zero()), MFe(Eigen::Array<double,
       3,1>::Zero()), MFd(Eigen::Array<double,3,1>::Zero()), MFu(Eigen::Array<
       double,3,1>::Zero()), MVWm(0), MVP(0), MVZ(0)

   , ZD(Eigen::Matrix<double,6,6>::Zero()), ZU(Eigen::Matrix<double,6,6>::Zero(
      )), ZE(Eigen::Matrix<double,6,6>::Zero()), ZV(Eigen::Matrix<double,6,6>
      ::Zero()), ZH(Eigen::Matrix<double,2,2>::Zero()), ZA(Eigen::Matrix<double,2,
      2>::Zero()), ZP(Eigen::Matrix<double,2,2>::Zero()), ZN(Eigen::Matrix<
      std::complex<double>,4,4>::Zero()), UV(Eigen::Matrix<std::complex<double>,6,
      6>::Zero()), UM(Eigen::Matrix<std::complex<double>,2,2>::Zero()), UP(
      Eigen::Matrix<std::complex<double>,2,2>::Zero()), ZEL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZER(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZDL(Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZDR(
      Eigen::Matrix<std::complex<double>,3,3>::Zero()), ZUL(Eigen::Matrix<
      std::complex<double>,3,3>::Zero()), ZUR(Eigen::Matrix<std::complex<double>,3
      ,3>::Zero()), ZZ(Eigen::Matrix<double,2,2>::Zero())

{
}

void MSSMRHN_physical::clear()
{
   MVG = 0.;
   MGlu = 0.;
   MSd = Eigen::Matrix<double,6,1>::Zero();
   ZD = Eigen::Matrix<double,6,6>::Zero();
   MSu = Eigen::Matrix<double,6,1>::Zero();
   ZU = Eigen::Matrix<double,6,6>::Zero();
   MSe = Eigen::Matrix<double,6,1>::Zero();
   ZE = Eigen::Matrix<double,6,6>::Zero();
   MSv = Eigen::Matrix<double,6,1>::Zero();
   ZV = Eigen::Matrix<double,6,6>::Zero();
   Mhh = Eigen::Matrix<double,2,1>::Zero();
   ZH = Eigen::Matrix<double,2,2>::Zero();
   MAh = Eigen::Matrix<double,2,1>::Zero();
   ZA = Eigen::Matrix<double,2,2>::Zero();
   MHpm = Eigen::Matrix<double,2,1>::Zero();
   ZP = Eigen::Matrix<double,2,2>::Zero();
   MChi = Eigen::Matrix<double,4,1>::Zero();
   ZN = Eigen::Matrix<std::complex<double>,4,4>::Zero();
   MFv = Eigen::Matrix<double,6,1>::Zero();
   UV = Eigen::Matrix<std::complex<double>,6,6>::Zero();
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

}

/**
 * Convert masses and mixing matrices to Haber-Kane convention:
 * Fermion masses are always positive and mixing matrices are allowed
 * to be complex.
 */
void MSSMRHN_physical::convert_to_hk()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));
   SLHA_io::convert_symmetric_fermion_mixings_to_hk(LOCALPHYSICAL(MFv), LOCALPHYSICAL(UV));

}

/**
 * Convert masses and mixing matrices to SLHA convention: Fermion
 * mixing matrices are always real and fermion masses are allowed to
 * be negative.
 */
void MSSMRHN_physical::convert_to_slha()
{
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MChi), LOCALPHYSICAL(ZN));
   SLHA_io::convert_symmetric_fermion_mixings_to_slha(LOCALPHYSICAL(MFv), LOCALPHYSICAL(UV));

}

Eigen::ArrayXd MSSMRHN_physical::get() const
{
   Eigen::ArrayXd pars(get_masses());

   pars.conservativeResize(444);

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
   pars(92) = ZU(0,0);
   pars(93) = ZU(0,1);
   pars(94) = ZU(0,2);
   pars(95) = ZU(0,3);
   pars(96) = ZU(0,4);
   pars(97) = ZU(0,5);
   pars(98) = ZU(1,0);
   pars(99) = ZU(1,1);
   pars(100) = ZU(1,2);
   pars(101) = ZU(1,3);
   pars(102) = ZU(1,4);
   pars(103) = ZU(1,5);
   pars(104) = ZU(2,0);
   pars(105) = ZU(2,1);
   pars(106) = ZU(2,2);
   pars(107) = ZU(2,3);
   pars(108) = ZU(2,4);
   pars(109) = ZU(2,5);
   pars(110) = ZU(3,0);
   pars(111) = ZU(3,1);
   pars(112) = ZU(3,2);
   pars(113) = ZU(3,3);
   pars(114) = ZU(3,4);
   pars(115) = ZU(3,5);
   pars(116) = ZU(4,0);
   pars(117) = ZU(4,1);
   pars(118) = ZU(4,2);
   pars(119) = ZU(4,3);
   pars(120) = ZU(4,4);
   pars(121) = ZU(4,5);
   pars(122) = ZU(5,0);
   pars(123) = ZU(5,1);
   pars(124) = ZU(5,2);
   pars(125) = ZU(5,3);
   pars(126) = ZU(5,4);
   pars(127) = ZU(5,5);
   pars(128) = ZE(0,0);
   pars(129) = ZE(0,1);
   pars(130) = ZE(0,2);
   pars(131) = ZE(0,3);
   pars(132) = ZE(0,4);
   pars(133) = ZE(0,5);
   pars(134) = ZE(1,0);
   pars(135) = ZE(1,1);
   pars(136) = ZE(1,2);
   pars(137) = ZE(1,3);
   pars(138) = ZE(1,4);
   pars(139) = ZE(1,5);
   pars(140) = ZE(2,0);
   pars(141) = ZE(2,1);
   pars(142) = ZE(2,2);
   pars(143) = ZE(2,3);
   pars(144) = ZE(2,4);
   pars(145) = ZE(2,5);
   pars(146) = ZE(3,0);
   pars(147) = ZE(3,1);
   pars(148) = ZE(3,2);
   pars(149) = ZE(3,3);
   pars(150) = ZE(3,4);
   pars(151) = ZE(3,5);
   pars(152) = ZE(4,0);
   pars(153) = ZE(4,1);
   pars(154) = ZE(4,2);
   pars(155) = ZE(4,3);
   pars(156) = ZE(4,4);
   pars(157) = ZE(4,5);
   pars(158) = ZE(5,0);
   pars(159) = ZE(5,1);
   pars(160) = ZE(5,2);
   pars(161) = ZE(5,3);
   pars(162) = ZE(5,4);
   pars(163) = ZE(5,5);
   pars(164) = ZV(0,0);
   pars(165) = ZV(0,1);
   pars(166) = ZV(0,2);
   pars(167) = ZV(0,3);
   pars(168) = ZV(0,4);
   pars(169) = ZV(0,5);
   pars(170) = ZV(1,0);
   pars(171) = ZV(1,1);
   pars(172) = ZV(1,2);
   pars(173) = ZV(1,3);
   pars(174) = ZV(1,4);
   pars(175) = ZV(1,5);
   pars(176) = ZV(2,0);
   pars(177) = ZV(2,1);
   pars(178) = ZV(2,2);
   pars(179) = ZV(2,3);
   pars(180) = ZV(2,4);
   pars(181) = ZV(2,5);
   pars(182) = ZV(3,0);
   pars(183) = ZV(3,1);
   pars(184) = ZV(3,2);
   pars(185) = ZV(3,3);
   pars(186) = ZV(3,4);
   pars(187) = ZV(3,5);
   pars(188) = ZV(4,0);
   pars(189) = ZV(4,1);
   pars(190) = ZV(4,2);
   pars(191) = ZV(4,3);
   pars(192) = ZV(4,4);
   pars(193) = ZV(4,5);
   pars(194) = ZV(5,0);
   pars(195) = ZV(5,1);
   pars(196) = ZV(5,2);
   pars(197) = ZV(5,3);
   pars(198) = ZV(5,4);
   pars(199) = ZV(5,5);
   pars(200) = ZH(0,0);
   pars(201) = ZH(0,1);
   pars(202) = ZH(1,0);
   pars(203) = ZH(1,1);
   pars(204) = ZA(0,0);
   pars(205) = ZA(0,1);
   pars(206) = ZA(1,0);
   pars(207) = ZA(1,1);
   pars(208) = ZP(0,0);
   pars(209) = ZP(0,1);
   pars(210) = ZP(1,0);
   pars(211) = ZP(1,1);
   pars(212) = Re(ZN(0,0));
   pars(213) = Im(ZN(0,0));
   pars(214) = Re(ZN(0,1));
   pars(215) = Im(ZN(0,1));
   pars(216) = Re(ZN(0,2));
   pars(217) = Im(ZN(0,2));
   pars(218) = Re(ZN(0,3));
   pars(219) = Im(ZN(0,3));
   pars(220) = Re(ZN(1,0));
   pars(221) = Im(ZN(1,0));
   pars(222) = Re(ZN(1,1));
   pars(223) = Im(ZN(1,1));
   pars(224) = Re(ZN(1,2));
   pars(225) = Im(ZN(1,2));
   pars(226) = Re(ZN(1,3));
   pars(227) = Im(ZN(1,3));
   pars(228) = Re(ZN(2,0));
   pars(229) = Im(ZN(2,0));
   pars(230) = Re(ZN(2,1));
   pars(231) = Im(ZN(2,1));
   pars(232) = Re(ZN(2,2));
   pars(233) = Im(ZN(2,2));
   pars(234) = Re(ZN(2,3));
   pars(235) = Im(ZN(2,3));
   pars(236) = Re(ZN(3,0));
   pars(237) = Im(ZN(3,0));
   pars(238) = Re(ZN(3,1));
   pars(239) = Im(ZN(3,1));
   pars(240) = Re(ZN(3,2));
   pars(241) = Im(ZN(3,2));
   pars(242) = Re(ZN(3,3));
   pars(243) = Im(ZN(3,3));
   pars(244) = Re(UV(0,0));
   pars(245) = Im(UV(0,0));
   pars(246) = Re(UV(0,1));
   pars(247) = Im(UV(0,1));
   pars(248) = Re(UV(0,2));
   pars(249) = Im(UV(0,2));
   pars(250) = Re(UV(0,3));
   pars(251) = Im(UV(0,3));
   pars(252) = Re(UV(0,4));
   pars(253) = Im(UV(0,4));
   pars(254) = Re(UV(0,5));
   pars(255) = Im(UV(0,5));
   pars(256) = Re(UV(1,0));
   pars(257) = Im(UV(1,0));
   pars(258) = Re(UV(1,1));
   pars(259) = Im(UV(1,1));
   pars(260) = Re(UV(1,2));
   pars(261) = Im(UV(1,2));
   pars(262) = Re(UV(1,3));
   pars(263) = Im(UV(1,3));
   pars(264) = Re(UV(1,4));
   pars(265) = Im(UV(1,4));
   pars(266) = Re(UV(1,5));
   pars(267) = Im(UV(1,5));
   pars(268) = Re(UV(2,0));
   pars(269) = Im(UV(2,0));
   pars(270) = Re(UV(2,1));
   pars(271) = Im(UV(2,1));
   pars(272) = Re(UV(2,2));
   pars(273) = Im(UV(2,2));
   pars(274) = Re(UV(2,3));
   pars(275) = Im(UV(2,3));
   pars(276) = Re(UV(2,4));
   pars(277) = Im(UV(2,4));
   pars(278) = Re(UV(2,5));
   pars(279) = Im(UV(2,5));
   pars(280) = Re(UV(3,0));
   pars(281) = Im(UV(3,0));
   pars(282) = Re(UV(3,1));
   pars(283) = Im(UV(3,1));
   pars(284) = Re(UV(3,2));
   pars(285) = Im(UV(3,2));
   pars(286) = Re(UV(3,3));
   pars(287) = Im(UV(3,3));
   pars(288) = Re(UV(3,4));
   pars(289) = Im(UV(3,4));
   pars(290) = Re(UV(3,5));
   pars(291) = Im(UV(3,5));
   pars(292) = Re(UV(4,0));
   pars(293) = Im(UV(4,0));
   pars(294) = Re(UV(4,1));
   pars(295) = Im(UV(4,1));
   pars(296) = Re(UV(4,2));
   pars(297) = Im(UV(4,2));
   pars(298) = Re(UV(4,3));
   pars(299) = Im(UV(4,3));
   pars(300) = Re(UV(4,4));
   pars(301) = Im(UV(4,4));
   pars(302) = Re(UV(4,5));
   pars(303) = Im(UV(4,5));
   pars(304) = Re(UV(5,0));
   pars(305) = Im(UV(5,0));
   pars(306) = Re(UV(5,1));
   pars(307) = Im(UV(5,1));
   pars(308) = Re(UV(5,2));
   pars(309) = Im(UV(5,2));
   pars(310) = Re(UV(5,3));
   pars(311) = Im(UV(5,3));
   pars(312) = Re(UV(5,4));
   pars(313) = Im(UV(5,4));
   pars(314) = Re(UV(5,5));
   pars(315) = Im(UV(5,5));
   pars(316) = Re(UM(0,0));
   pars(317) = Im(UM(0,0));
   pars(318) = Re(UM(0,1));
   pars(319) = Im(UM(0,1));
   pars(320) = Re(UM(1,0));
   pars(321) = Im(UM(1,0));
   pars(322) = Re(UM(1,1));
   pars(323) = Im(UM(1,1));
   pars(324) = Re(UP(0,0));
   pars(325) = Im(UP(0,0));
   pars(326) = Re(UP(0,1));
   pars(327) = Im(UP(0,1));
   pars(328) = Re(UP(1,0));
   pars(329) = Im(UP(1,0));
   pars(330) = Re(UP(1,1));
   pars(331) = Im(UP(1,1));
   pars(332) = Re(ZEL(0,0));
   pars(333) = Im(ZEL(0,0));
   pars(334) = Re(ZEL(0,1));
   pars(335) = Im(ZEL(0,1));
   pars(336) = Re(ZEL(0,2));
   pars(337) = Im(ZEL(0,2));
   pars(338) = Re(ZEL(1,0));
   pars(339) = Im(ZEL(1,0));
   pars(340) = Re(ZEL(1,1));
   pars(341) = Im(ZEL(1,1));
   pars(342) = Re(ZEL(1,2));
   pars(343) = Im(ZEL(1,2));
   pars(344) = Re(ZEL(2,0));
   pars(345) = Im(ZEL(2,0));
   pars(346) = Re(ZEL(2,1));
   pars(347) = Im(ZEL(2,1));
   pars(348) = Re(ZEL(2,2));
   pars(349) = Im(ZEL(2,2));
   pars(350) = Re(ZER(0,0));
   pars(351) = Im(ZER(0,0));
   pars(352) = Re(ZER(0,1));
   pars(353) = Im(ZER(0,1));
   pars(354) = Re(ZER(0,2));
   pars(355) = Im(ZER(0,2));
   pars(356) = Re(ZER(1,0));
   pars(357) = Im(ZER(1,0));
   pars(358) = Re(ZER(1,1));
   pars(359) = Im(ZER(1,1));
   pars(360) = Re(ZER(1,2));
   pars(361) = Im(ZER(1,2));
   pars(362) = Re(ZER(2,0));
   pars(363) = Im(ZER(2,0));
   pars(364) = Re(ZER(2,1));
   pars(365) = Im(ZER(2,1));
   pars(366) = Re(ZER(2,2));
   pars(367) = Im(ZER(2,2));
   pars(368) = Re(ZDL(0,0));
   pars(369) = Im(ZDL(0,0));
   pars(370) = Re(ZDL(0,1));
   pars(371) = Im(ZDL(0,1));
   pars(372) = Re(ZDL(0,2));
   pars(373) = Im(ZDL(0,2));
   pars(374) = Re(ZDL(1,0));
   pars(375) = Im(ZDL(1,0));
   pars(376) = Re(ZDL(1,1));
   pars(377) = Im(ZDL(1,1));
   pars(378) = Re(ZDL(1,2));
   pars(379) = Im(ZDL(1,2));
   pars(380) = Re(ZDL(2,0));
   pars(381) = Im(ZDL(2,0));
   pars(382) = Re(ZDL(2,1));
   pars(383) = Im(ZDL(2,1));
   pars(384) = Re(ZDL(2,2));
   pars(385) = Im(ZDL(2,2));
   pars(386) = Re(ZDR(0,0));
   pars(387) = Im(ZDR(0,0));
   pars(388) = Re(ZDR(0,1));
   pars(389) = Im(ZDR(0,1));
   pars(390) = Re(ZDR(0,2));
   pars(391) = Im(ZDR(0,2));
   pars(392) = Re(ZDR(1,0));
   pars(393) = Im(ZDR(1,0));
   pars(394) = Re(ZDR(1,1));
   pars(395) = Im(ZDR(1,1));
   pars(396) = Re(ZDR(1,2));
   pars(397) = Im(ZDR(1,2));
   pars(398) = Re(ZDR(2,0));
   pars(399) = Im(ZDR(2,0));
   pars(400) = Re(ZDR(2,1));
   pars(401) = Im(ZDR(2,1));
   pars(402) = Re(ZDR(2,2));
   pars(403) = Im(ZDR(2,2));
   pars(404) = Re(ZUL(0,0));
   pars(405) = Im(ZUL(0,0));
   pars(406) = Re(ZUL(0,1));
   pars(407) = Im(ZUL(0,1));
   pars(408) = Re(ZUL(0,2));
   pars(409) = Im(ZUL(0,2));
   pars(410) = Re(ZUL(1,0));
   pars(411) = Im(ZUL(1,0));
   pars(412) = Re(ZUL(1,1));
   pars(413) = Im(ZUL(1,1));
   pars(414) = Re(ZUL(1,2));
   pars(415) = Im(ZUL(1,2));
   pars(416) = Re(ZUL(2,0));
   pars(417) = Im(ZUL(2,0));
   pars(418) = Re(ZUL(2,1));
   pars(419) = Im(ZUL(2,1));
   pars(420) = Re(ZUL(2,2));
   pars(421) = Im(ZUL(2,2));
   pars(422) = Re(ZUR(0,0));
   pars(423) = Im(ZUR(0,0));
   pars(424) = Re(ZUR(0,1));
   pars(425) = Im(ZUR(0,1));
   pars(426) = Re(ZUR(0,2));
   pars(427) = Im(ZUR(0,2));
   pars(428) = Re(ZUR(1,0));
   pars(429) = Im(ZUR(1,0));
   pars(430) = Re(ZUR(1,1));
   pars(431) = Im(ZUR(1,1));
   pars(432) = Re(ZUR(1,2));
   pars(433) = Im(ZUR(1,2));
   pars(434) = Re(ZUR(2,0));
   pars(435) = Im(ZUR(2,0));
   pars(436) = Re(ZUR(2,1));
   pars(437) = Im(ZUR(2,1));
   pars(438) = Re(ZUR(2,2));
   pars(439) = Im(ZUR(2,2));
   pars(440) = ZZ(0,0);
   pars(441) = ZZ(0,1);
   pars(442) = ZZ(1,0);
   pars(443) = ZZ(1,1);


   return pars;
}

void MSSMRHN_physical::set(const Eigen::ArrayXd& pars)
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
   ZU(0,0) = pars(92);
   ZU(0,1) = pars(93);
   ZU(0,2) = pars(94);
   ZU(0,3) = pars(95);
   ZU(0,4) = pars(96);
   ZU(0,5) = pars(97);
   ZU(1,0) = pars(98);
   ZU(1,1) = pars(99);
   ZU(1,2) = pars(100);
   ZU(1,3) = pars(101);
   ZU(1,4) = pars(102);
   ZU(1,5) = pars(103);
   ZU(2,0) = pars(104);
   ZU(2,1) = pars(105);
   ZU(2,2) = pars(106);
   ZU(2,3) = pars(107);
   ZU(2,4) = pars(108);
   ZU(2,5) = pars(109);
   ZU(3,0) = pars(110);
   ZU(3,1) = pars(111);
   ZU(3,2) = pars(112);
   ZU(3,3) = pars(113);
   ZU(3,4) = pars(114);
   ZU(3,5) = pars(115);
   ZU(4,0) = pars(116);
   ZU(4,1) = pars(117);
   ZU(4,2) = pars(118);
   ZU(4,3) = pars(119);
   ZU(4,4) = pars(120);
   ZU(4,5) = pars(121);
   ZU(5,0) = pars(122);
   ZU(5,1) = pars(123);
   ZU(5,2) = pars(124);
   ZU(5,3) = pars(125);
   ZU(5,4) = pars(126);
   ZU(5,5) = pars(127);
   ZE(0,0) = pars(128);
   ZE(0,1) = pars(129);
   ZE(0,2) = pars(130);
   ZE(0,3) = pars(131);
   ZE(0,4) = pars(132);
   ZE(0,5) = pars(133);
   ZE(1,0) = pars(134);
   ZE(1,1) = pars(135);
   ZE(1,2) = pars(136);
   ZE(1,3) = pars(137);
   ZE(1,4) = pars(138);
   ZE(1,5) = pars(139);
   ZE(2,0) = pars(140);
   ZE(2,1) = pars(141);
   ZE(2,2) = pars(142);
   ZE(2,3) = pars(143);
   ZE(2,4) = pars(144);
   ZE(2,5) = pars(145);
   ZE(3,0) = pars(146);
   ZE(3,1) = pars(147);
   ZE(3,2) = pars(148);
   ZE(3,3) = pars(149);
   ZE(3,4) = pars(150);
   ZE(3,5) = pars(151);
   ZE(4,0) = pars(152);
   ZE(4,1) = pars(153);
   ZE(4,2) = pars(154);
   ZE(4,3) = pars(155);
   ZE(4,4) = pars(156);
   ZE(4,5) = pars(157);
   ZE(5,0) = pars(158);
   ZE(5,1) = pars(159);
   ZE(5,2) = pars(160);
   ZE(5,3) = pars(161);
   ZE(5,4) = pars(162);
   ZE(5,5) = pars(163);
   ZV(0,0) = pars(164);
   ZV(0,1) = pars(165);
   ZV(0,2) = pars(166);
   ZV(0,3) = pars(167);
   ZV(0,4) = pars(168);
   ZV(0,5) = pars(169);
   ZV(1,0) = pars(170);
   ZV(1,1) = pars(171);
   ZV(1,2) = pars(172);
   ZV(1,3) = pars(173);
   ZV(1,4) = pars(174);
   ZV(1,5) = pars(175);
   ZV(2,0) = pars(176);
   ZV(2,1) = pars(177);
   ZV(2,2) = pars(178);
   ZV(2,3) = pars(179);
   ZV(2,4) = pars(180);
   ZV(2,5) = pars(181);
   ZV(3,0) = pars(182);
   ZV(3,1) = pars(183);
   ZV(3,2) = pars(184);
   ZV(3,3) = pars(185);
   ZV(3,4) = pars(186);
   ZV(3,5) = pars(187);
   ZV(4,0) = pars(188);
   ZV(4,1) = pars(189);
   ZV(4,2) = pars(190);
   ZV(4,3) = pars(191);
   ZV(4,4) = pars(192);
   ZV(4,5) = pars(193);
   ZV(5,0) = pars(194);
   ZV(5,1) = pars(195);
   ZV(5,2) = pars(196);
   ZV(5,3) = pars(197);
   ZV(5,4) = pars(198);
   ZV(5,5) = pars(199);
   ZH(0,0) = pars(200);
   ZH(0,1) = pars(201);
   ZH(1,0) = pars(202);
   ZH(1,1) = pars(203);
   ZA(0,0) = pars(204);
   ZA(0,1) = pars(205);
   ZA(1,0) = pars(206);
   ZA(1,1) = pars(207);
   ZP(0,0) = pars(208);
   ZP(0,1) = pars(209);
   ZP(1,0) = pars(210);
   ZP(1,1) = pars(211);
   ZN(0,0) = std::complex<double>(pars(212), pars(213));
   ZN(0,1) = std::complex<double>(pars(214), pars(215));
   ZN(0,2) = std::complex<double>(pars(216), pars(217));
   ZN(0,3) = std::complex<double>(pars(218), pars(219));
   ZN(1,0) = std::complex<double>(pars(220), pars(221));
   ZN(1,1) = std::complex<double>(pars(222), pars(223));
   ZN(1,2) = std::complex<double>(pars(224), pars(225));
   ZN(1,3) = std::complex<double>(pars(226), pars(227));
   ZN(2,0) = std::complex<double>(pars(228), pars(229));
   ZN(2,1) = std::complex<double>(pars(230), pars(231));
   ZN(2,2) = std::complex<double>(pars(232), pars(233));
   ZN(2,3) = std::complex<double>(pars(234), pars(235));
   ZN(3,0) = std::complex<double>(pars(236), pars(237));
   ZN(3,1) = std::complex<double>(pars(238), pars(239));
   ZN(3,2) = std::complex<double>(pars(240), pars(241));
   ZN(3,3) = std::complex<double>(pars(242), pars(243));
   UV(0,0) = std::complex<double>(pars(244), pars(245));
   UV(0,1) = std::complex<double>(pars(246), pars(247));
   UV(0,2) = std::complex<double>(pars(248), pars(249));
   UV(0,3) = std::complex<double>(pars(250), pars(251));
   UV(0,4) = std::complex<double>(pars(252), pars(253));
   UV(0,5) = std::complex<double>(pars(254), pars(255));
   UV(1,0) = std::complex<double>(pars(256), pars(257));
   UV(1,1) = std::complex<double>(pars(258), pars(259));
   UV(1,2) = std::complex<double>(pars(260), pars(261));
   UV(1,3) = std::complex<double>(pars(262), pars(263));
   UV(1,4) = std::complex<double>(pars(264), pars(265));
   UV(1,5) = std::complex<double>(pars(266), pars(267));
   UV(2,0) = std::complex<double>(pars(268), pars(269));
   UV(2,1) = std::complex<double>(pars(270), pars(271));
   UV(2,2) = std::complex<double>(pars(272), pars(273));
   UV(2,3) = std::complex<double>(pars(274), pars(275));
   UV(2,4) = std::complex<double>(pars(276), pars(277));
   UV(2,5) = std::complex<double>(pars(278), pars(279));
   UV(3,0) = std::complex<double>(pars(280), pars(281));
   UV(3,1) = std::complex<double>(pars(282), pars(283));
   UV(3,2) = std::complex<double>(pars(284), pars(285));
   UV(3,3) = std::complex<double>(pars(286), pars(287));
   UV(3,4) = std::complex<double>(pars(288), pars(289));
   UV(3,5) = std::complex<double>(pars(290), pars(291));
   UV(4,0) = std::complex<double>(pars(292), pars(293));
   UV(4,1) = std::complex<double>(pars(294), pars(295));
   UV(4,2) = std::complex<double>(pars(296), pars(297));
   UV(4,3) = std::complex<double>(pars(298), pars(299));
   UV(4,4) = std::complex<double>(pars(300), pars(301));
   UV(4,5) = std::complex<double>(pars(302), pars(303));
   UV(5,0) = std::complex<double>(pars(304), pars(305));
   UV(5,1) = std::complex<double>(pars(306), pars(307));
   UV(5,2) = std::complex<double>(pars(308), pars(309));
   UV(5,3) = std::complex<double>(pars(310), pars(311));
   UV(5,4) = std::complex<double>(pars(312), pars(313));
   UV(5,5) = std::complex<double>(pars(314), pars(315));
   UM(0,0) = std::complex<double>(pars(316), pars(317));
   UM(0,1) = std::complex<double>(pars(318), pars(319));
   UM(1,0) = std::complex<double>(pars(320), pars(321));
   UM(1,1) = std::complex<double>(pars(322), pars(323));
   UP(0,0) = std::complex<double>(pars(324), pars(325));
   UP(0,1) = std::complex<double>(pars(326), pars(327));
   UP(1,0) = std::complex<double>(pars(328), pars(329));
   UP(1,1) = std::complex<double>(pars(330), pars(331));
   ZEL(0,0) = std::complex<double>(pars(332), pars(333));
   ZEL(0,1) = std::complex<double>(pars(334), pars(335));
   ZEL(0,2) = std::complex<double>(pars(336), pars(337));
   ZEL(1,0) = std::complex<double>(pars(338), pars(339));
   ZEL(1,1) = std::complex<double>(pars(340), pars(341));
   ZEL(1,2) = std::complex<double>(pars(342), pars(343));
   ZEL(2,0) = std::complex<double>(pars(344), pars(345));
   ZEL(2,1) = std::complex<double>(pars(346), pars(347));
   ZEL(2,2) = std::complex<double>(pars(348), pars(349));
   ZER(0,0) = std::complex<double>(pars(350), pars(351));
   ZER(0,1) = std::complex<double>(pars(352), pars(353));
   ZER(0,2) = std::complex<double>(pars(354), pars(355));
   ZER(1,0) = std::complex<double>(pars(356), pars(357));
   ZER(1,1) = std::complex<double>(pars(358), pars(359));
   ZER(1,2) = std::complex<double>(pars(360), pars(361));
   ZER(2,0) = std::complex<double>(pars(362), pars(363));
   ZER(2,1) = std::complex<double>(pars(364), pars(365));
   ZER(2,2) = std::complex<double>(pars(366), pars(367));
   ZDL(0,0) = std::complex<double>(pars(368), pars(369));
   ZDL(0,1) = std::complex<double>(pars(370), pars(371));
   ZDL(0,2) = std::complex<double>(pars(372), pars(373));
   ZDL(1,0) = std::complex<double>(pars(374), pars(375));
   ZDL(1,1) = std::complex<double>(pars(376), pars(377));
   ZDL(1,2) = std::complex<double>(pars(378), pars(379));
   ZDL(2,0) = std::complex<double>(pars(380), pars(381));
   ZDL(2,1) = std::complex<double>(pars(382), pars(383));
   ZDL(2,2) = std::complex<double>(pars(384), pars(385));
   ZDR(0,0) = std::complex<double>(pars(386), pars(387));
   ZDR(0,1) = std::complex<double>(pars(388), pars(389));
   ZDR(0,2) = std::complex<double>(pars(390), pars(391));
   ZDR(1,0) = std::complex<double>(pars(392), pars(393));
   ZDR(1,1) = std::complex<double>(pars(394), pars(395));
   ZDR(1,2) = std::complex<double>(pars(396), pars(397));
   ZDR(2,0) = std::complex<double>(pars(398), pars(399));
   ZDR(2,1) = std::complex<double>(pars(400), pars(401));
   ZDR(2,2) = std::complex<double>(pars(402), pars(403));
   ZUL(0,0) = std::complex<double>(pars(404), pars(405));
   ZUL(0,1) = std::complex<double>(pars(406), pars(407));
   ZUL(0,2) = std::complex<double>(pars(408), pars(409));
   ZUL(1,0) = std::complex<double>(pars(410), pars(411));
   ZUL(1,1) = std::complex<double>(pars(412), pars(413));
   ZUL(1,2) = std::complex<double>(pars(414), pars(415));
   ZUL(2,0) = std::complex<double>(pars(416), pars(417));
   ZUL(2,1) = std::complex<double>(pars(418), pars(419));
   ZUL(2,2) = std::complex<double>(pars(420), pars(421));
   ZUR(0,0) = std::complex<double>(pars(422), pars(423));
   ZUR(0,1) = std::complex<double>(pars(424), pars(425));
   ZUR(0,2) = std::complex<double>(pars(426), pars(427));
   ZUR(1,0) = std::complex<double>(pars(428), pars(429));
   ZUR(1,1) = std::complex<double>(pars(430), pars(431));
   ZUR(1,2) = std::complex<double>(pars(432), pars(433));
   ZUR(2,0) = std::complex<double>(pars(434), pars(435));
   ZUR(2,1) = std::complex<double>(pars(436), pars(437));
   ZUR(2,2) = std::complex<double>(pars(438), pars(439));
   ZZ(0,0) = pars(440);
   ZZ(0,1) = pars(441);
   ZZ(1,0) = pars(442);
   ZZ(1,1) = pars(443);

}

Eigen::ArrayXd MSSMRHN_physical::get_masses() const
{
   Eigen::ArrayXd pars(56);

   pars(0) = MVG;
   pars(1) = MGlu;
   pars(2) = MSd(0);
   pars(3) = MSd(1);
   pars(4) = MSd(2);
   pars(5) = MSd(3);
   pars(6) = MSd(4);
   pars(7) = MSd(5);
   pars(8) = MSu(0);
   pars(9) = MSu(1);
   pars(10) = MSu(2);
   pars(11) = MSu(3);
   pars(12) = MSu(4);
   pars(13) = MSu(5);
   pars(14) = MSe(0);
   pars(15) = MSe(1);
   pars(16) = MSe(2);
   pars(17) = MSe(3);
   pars(18) = MSe(4);
   pars(19) = MSe(5);
   pars(20) = MSv(0);
   pars(21) = MSv(1);
   pars(22) = MSv(2);
   pars(23) = MSv(3);
   pars(24) = MSv(4);
   pars(25) = MSv(5);
   pars(26) = Mhh(0);
   pars(27) = Mhh(1);
   pars(28) = MAh(0);
   pars(29) = MAh(1);
   pars(30) = MHpm(0);
   pars(31) = MHpm(1);
   pars(32) = MChi(0);
   pars(33) = MChi(1);
   pars(34) = MChi(2);
   pars(35) = MChi(3);
   pars(36) = MFv(0);
   pars(37) = MFv(1);
   pars(38) = MFv(2);
   pars(39) = MFv(3);
   pars(40) = MFv(4);
   pars(41) = MFv(5);
   pars(42) = MCha(0);
   pars(43) = MCha(1);
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

void MSSMRHN_physical::set_masses(const Eigen::ArrayXd& pars)
{
   MVG = pars(0);
   MGlu = pars(1);
   MSd(0) = pars(2);
   MSd(1) = pars(3);
   MSd(2) = pars(4);
   MSd(3) = pars(5);
   MSd(4) = pars(6);
   MSd(5) = pars(7);
   MSu(0) = pars(8);
   MSu(1) = pars(9);
   MSu(2) = pars(10);
   MSu(3) = pars(11);
   MSu(4) = pars(12);
   MSu(5) = pars(13);
   MSe(0) = pars(14);
   MSe(1) = pars(15);
   MSe(2) = pars(16);
   MSe(3) = pars(17);
   MSe(4) = pars(18);
   MSe(5) = pars(19);
   MSv(0) = pars(20);
   MSv(1) = pars(21);
   MSv(2) = pars(22);
   MSv(3) = pars(23);
   MSv(4) = pars(24);
   MSv(5) = pars(25);
   Mhh(0) = pars(26);
   Mhh(1) = pars(27);
   MAh(0) = pars(28);
   MAh(1) = pars(29);
   MHpm(0) = pars(30);
   MHpm(1) = pars(31);
   MChi(0) = pars(32);
   MChi(1) = pars(33);
   MChi(2) = pars(34);
   MChi(3) = pars(35);
   MFv(0) = pars(36);
   MFv(1) = pars(37);
   MFv(2) = pars(38);
   MFv(3) = pars(39);
   MFv(4) = pars(40);
   MFv(5) = pars(41);
   MCha(0) = pars(42);
   MCha(1) = pars(43);
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

void MSSMRHN_physical::print(std::ostream& ostr) const
{
   ostr << "----------------------------------------\n"
           "pole masses:\n"
           "----------------------------------------\n";
   ostr << "MVG = " << MVG << '\n';
   ostr << "MGlu = " << MGlu << '\n';
   ostr << "MSd = " << MSd.transpose() << '\n';
   ostr << "MSu = " << MSu.transpose() << '\n';
   ostr << "MSe = " << MSe.transpose() << '\n';
   ostr << "MSv = " << MSv.transpose() << '\n';
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
   ostr << "ZU = " << ZU << '\n';
   ostr << "ZE = " << ZE << '\n';
   ostr << "ZV = " << ZV << '\n';
   ostr << "ZH = " << ZH << '\n';
   ostr << "ZA = " << ZA << '\n';
   ostr << "ZP = " << ZP << '\n';
   ostr << "ZN = " << ZN << '\n';
   ostr << "UV = " << UV << '\n';
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

std::ostream& operator<<(std::ostream& ostr, const MSSMRHN_physical& phys_pars)
{
   phys_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
