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


#include "E6SSM_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES(l) calc_soft_traces(l);

const int E6SSM_soft_parameters::numberOfParameters;

E6SSM_soft_parameters::E6SSM_soft_parameters(const E6SSM_input_parameters& input_)
   : E6SSM_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

E6SSM_soft_parameters::E6SSM_soft_parameters(
   const E6SSM_susy_parameters& susy_model
   , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>& TYe_,
   const Eigen::Matrix<double,3,3>& TKappa_, const Eigen::Matrix<double,2,2>&
   TLambda12_, double TLambdax_, const Eigen::Matrix<double,3,3>& TYu_, double
   BMuPr_, const Eigen::Matrix<double,3,3>& mq2_, const Eigen::Matrix<double,3,
   3>& ml2_, double mHd2_, double mHu2_, const Eigen::Matrix<double,3,3>& md2_,
   const Eigen::Matrix<double,3,3>& mu2_, const Eigen::Matrix<double,3,3>& me2_
   , double ms2_, const Eigen::Matrix<double,2,2>& mH1I2_, const Eigen::Matrix<
   double,2,2>& mH2I2_, const Eigen::Matrix<double,2,2>& msI2_, const Eigen::
   Matrix<double,3,3>& mDx2_, const Eigen::Matrix<double,3,3>& mDxbar2_, double
    mHp2_, double mHpbar2_, double MassB_, double MassWB_, double MassG_,
   double MassBp_
)
   : E6SSM_susy_parameters(susy_model)
   , TYd(TYd_), TYe(TYe_), TKappa(TKappa_), TLambda12(TLambda12_), TLambdax(
   TLambdax_), TYu(TYu_), BMuPr(BMuPr_), mq2(mq2_), ml2(ml2_), mHd2(mHd2_),
   mHu2(mHu2_), md2(md2_), mu2(mu2_), me2(me2_), ms2(ms2_), mH1I2(mH1I2_),
   mH2I2(mH2I2_), msI2(msI2_), mDx2(mDx2_), mDxbar2(mDxbar2_), mHp2(mHp2_),
   mHpbar2(mHpbar2_), MassB(MassB_), MassWB(MassWB_), MassG(MassG_), MassBp(
   MassBp_)
{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd E6SSM_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

E6SSM_soft_parameters E6SSM_soft_parameters::calc_beta(int loops) const
{
   Eigen::Matrix<double,3,3> beta_TYd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_TYe = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_TKappa = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,2,2> beta_TLambda12 = Eigen::Matrix<double,2,2>::Zero();
   double beta_TLambdax = 0.;
   Eigen::Matrix<double,3,3> beta_TYu = Eigen::Matrix<double,3,3>::Zero();
   double beta_BMuPr = 0.;
   Eigen::Matrix<double,3,3> beta_mq2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_ml2 = Eigen::Matrix<double,3,3>::Zero();
   double beta_mHd2 = 0.;
   double beta_mHu2 = 0.;
   Eigen::Matrix<double,3,3> beta_md2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_mu2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_me2 = Eigen::Matrix<double,3,3>::Zero();
   double beta_ms2 = 0.;
   Eigen::Matrix<double,2,2> beta_mH1I2 = Eigen::Matrix<double,2,2>::Zero();
   Eigen::Matrix<double,2,2> beta_mH2I2 = Eigen::Matrix<double,2,2>::Zero();
   Eigen::Matrix<double,2,2> beta_msI2 = Eigen::Matrix<double,2,2>::Zero();
   Eigen::Matrix<double,3,3> beta_mDx2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_mDxbar2 = Eigen::Matrix<double,3,3>::Zero();
   double beta_mHp2 = 0.;
   double beta_mHpbar2 = 0.;
   double beta_MassB = 0.;
   double beta_MassWB = 0.;
   double beta_MassG = 0.;
   double beta_MassBp = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_TYd += calc_beta_TYd_1_loop(TRACE_STRUCT);
      beta_TYe += calc_beta_TYe_1_loop(TRACE_STRUCT);
      beta_TKappa += calc_beta_TKappa_1_loop(TRACE_STRUCT);
      beta_TLambda12 += calc_beta_TLambda12_1_loop(TRACE_STRUCT);
      beta_TLambdax += calc_beta_TLambdax_1_loop(TRACE_STRUCT);
      beta_TYu += calc_beta_TYu_1_loop(TRACE_STRUCT);
      beta_BMuPr += calc_beta_BMuPr_1_loop(TRACE_STRUCT);
      beta_mq2 += calc_beta_mq2_1_loop(TRACE_STRUCT);
      beta_ml2 += calc_beta_ml2_1_loop(TRACE_STRUCT);
      beta_mHd2 += calc_beta_mHd2_1_loop(TRACE_STRUCT);
      beta_mHu2 += calc_beta_mHu2_1_loop(TRACE_STRUCT);
      beta_md2 += calc_beta_md2_1_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_1_loop(TRACE_STRUCT);
      beta_me2 += calc_beta_me2_1_loop(TRACE_STRUCT);
      beta_ms2 += calc_beta_ms2_1_loop(TRACE_STRUCT);
      beta_mH1I2 += calc_beta_mH1I2_1_loop(TRACE_STRUCT);
      beta_mH2I2 += calc_beta_mH2I2_1_loop(TRACE_STRUCT);
      beta_msI2 += calc_beta_msI2_1_loop(TRACE_STRUCT);
      beta_mDx2 += calc_beta_mDx2_1_loop(TRACE_STRUCT);
      beta_mDxbar2 += calc_beta_mDxbar2_1_loop(TRACE_STRUCT);
      beta_mHp2 += calc_beta_mHp2_1_loop(TRACE_STRUCT);
      beta_mHpbar2 += calc_beta_mHpbar2_1_loop(TRACE_STRUCT);
      beta_MassB += calc_beta_MassB_1_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_1_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_1_loop(TRACE_STRUCT);
      beta_MassBp += calc_beta_MassBp_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_TYd += calc_beta_TYd_2_loop(TRACE_STRUCT);
         beta_TYe += calc_beta_TYe_2_loop(TRACE_STRUCT);
         beta_TKappa += calc_beta_TKappa_2_loop(TRACE_STRUCT);
         beta_TLambda12 += calc_beta_TLambda12_2_loop(TRACE_STRUCT);
         beta_TLambdax += calc_beta_TLambdax_2_loop(TRACE_STRUCT);
         beta_TYu += calc_beta_TYu_2_loop(TRACE_STRUCT);
         beta_BMuPr += calc_beta_BMuPr_2_loop(TRACE_STRUCT);
         beta_mq2 += calc_beta_mq2_2_loop(TRACE_STRUCT);
         beta_ml2 += calc_beta_ml2_2_loop(TRACE_STRUCT);
         beta_mHd2 += calc_beta_mHd2_2_loop(TRACE_STRUCT);
         beta_mHu2 += calc_beta_mHu2_2_loop(TRACE_STRUCT);
         beta_md2 += calc_beta_md2_2_loop(TRACE_STRUCT);
         beta_mu2 += calc_beta_mu2_2_loop(TRACE_STRUCT);
         beta_me2 += calc_beta_me2_2_loop(TRACE_STRUCT);
         beta_ms2 += calc_beta_ms2_2_loop(TRACE_STRUCT);
         beta_mH1I2 += calc_beta_mH1I2_2_loop(TRACE_STRUCT);
         beta_mH2I2 += calc_beta_mH2I2_2_loop(TRACE_STRUCT);
         beta_msI2 += calc_beta_msI2_2_loop(TRACE_STRUCT);
         beta_mDx2 += calc_beta_mDx2_2_loop(TRACE_STRUCT);
         beta_mDxbar2 += calc_beta_mDxbar2_2_loop(TRACE_STRUCT);
         beta_mHp2 += calc_beta_mHp2_2_loop(TRACE_STRUCT);
         beta_mHpbar2 += calc_beta_mHpbar2_2_loop(TRACE_STRUCT);
         beta_MassB += calc_beta_MassB_2_loop(TRACE_STRUCT);
         beta_MassWB += calc_beta_MassWB_2_loop(TRACE_STRUCT);
         beta_MassG += calc_beta_MassG_2_loop(TRACE_STRUCT);
         beta_MassBp += calc_beta_MassBp_2_loop(TRACE_STRUCT);

         if (loops > 2) {

            if (loops > 3) {

               if (loops > 4) {

               }
            }
         }
      }
   }


   const E6SSM_susy_parameters susy_betas(E6SSM_susy_parameters::calc_beta(loops));

   return E6SSM_soft_parameters(susy_betas, beta_TYd, beta_TYe, beta_TKappa, beta_TLambda12, beta_TLambdax, beta_TYu, beta_BMuPr, beta_mq2, beta_ml2, beta_mHd2, beta_mHu2, beta_md2, beta_mu2, beta_me2, beta_ms2, beta_mH1I2, beta_mH2I2, beta_msI2, beta_mDx2, beta_mDxbar2, beta_mHp2, beta_mHpbar2, beta_MassB, beta_MassWB, beta_MassG, beta_MassBp);
}

E6SSM_soft_parameters E6SSM_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void E6SSM_soft_parameters::clear()
{
   E6SSM_susy_parameters::clear();

   TYd = Eigen::Matrix<double,3,3>::Zero();
   TYe = Eigen::Matrix<double,3,3>::Zero();
   TKappa = Eigen::Matrix<double,3,3>::Zero();
   TLambda12 = Eigen::Matrix<double,2,2>::Zero();
   TLambdax = 0.;
   TYu = Eigen::Matrix<double,3,3>::Zero();
   BMuPr = 0.;
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   ms2 = 0.;
   mH1I2 = Eigen::Matrix<double,2,2>::Zero();
   mH2I2 = Eigen::Matrix<double,2,2>::Zero();
   msI2 = Eigen::Matrix<double,2,2>::Zero();
   mDx2 = Eigen::Matrix<double,3,3>::Zero();
   mDxbar2 = Eigen::Matrix<double,3,3>::Zero();
   mHp2 = 0.;
   mHpbar2 = 0.;
   MassB = 0.;
   MassWB = 0.;
   MassG = 0.;
   MassBp = 0.;

}

Eigen::ArrayXd E6SSM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(E6SSM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(49) = TYd(0,0);
   pars(50) = TYd(0,1);
   pars(51) = TYd(0,2);
   pars(52) = TYd(1,0);
   pars(53) = TYd(1,1);
   pars(54) = TYd(1,2);
   pars(55) = TYd(2,0);
   pars(56) = TYd(2,1);
   pars(57) = TYd(2,2);
   pars(58) = TYe(0,0);
   pars(59) = TYe(0,1);
   pars(60) = TYe(0,2);
   pars(61) = TYe(1,0);
   pars(62) = TYe(1,1);
   pars(63) = TYe(1,2);
   pars(64) = TYe(2,0);
   pars(65) = TYe(2,1);
   pars(66) = TYe(2,2);
   pars(67) = TKappa(0,0);
   pars(68) = TKappa(0,1);
   pars(69) = TKappa(0,2);
   pars(70) = TKappa(1,0);
   pars(71) = TKappa(1,1);
   pars(72) = TKappa(1,2);
   pars(73) = TKappa(2,0);
   pars(74) = TKappa(2,1);
   pars(75) = TKappa(2,2);
   pars(76) = TLambda12(0,0);
   pars(77) = TLambda12(0,1);
   pars(78) = TLambda12(1,0);
   pars(79) = TLambda12(1,1);
   pars(80) = TLambdax;
   pars(81) = TYu(0,0);
   pars(82) = TYu(0,1);
   pars(83) = TYu(0,2);
   pars(84) = TYu(1,0);
   pars(85) = TYu(1,1);
   pars(86) = TYu(1,2);
   pars(87) = TYu(2,0);
   pars(88) = TYu(2,1);
   pars(89) = TYu(2,2);
   pars(90) = BMuPr;
   pars(91) = mq2(0,0);
   pars(92) = mq2(0,1);
   pars(93) = mq2(0,2);
   pars(94) = mq2(1,0);
   pars(95) = mq2(1,1);
   pars(96) = mq2(1,2);
   pars(97) = mq2(2,0);
   pars(98) = mq2(2,1);
   pars(99) = mq2(2,2);
   pars(100) = ml2(0,0);
   pars(101) = ml2(0,1);
   pars(102) = ml2(0,2);
   pars(103) = ml2(1,0);
   pars(104) = ml2(1,1);
   pars(105) = ml2(1,2);
   pars(106) = ml2(2,0);
   pars(107) = ml2(2,1);
   pars(108) = ml2(2,2);
   pars(109) = mHd2;
   pars(110) = mHu2;
   pars(111) = md2(0,0);
   pars(112) = md2(0,1);
   pars(113) = md2(0,2);
   pars(114) = md2(1,0);
   pars(115) = md2(1,1);
   pars(116) = md2(1,2);
   pars(117) = md2(2,0);
   pars(118) = md2(2,1);
   pars(119) = md2(2,2);
   pars(120) = mu2(0,0);
   pars(121) = mu2(0,1);
   pars(122) = mu2(0,2);
   pars(123) = mu2(1,0);
   pars(124) = mu2(1,1);
   pars(125) = mu2(1,2);
   pars(126) = mu2(2,0);
   pars(127) = mu2(2,1);
   pars(128) = mu2(2,2);
   pars(129) = me2(0,0);
   pars(130) = me2(0,1);
   pars(131) = me2(0,2);
   pars(132) = me2(1,0);
   pars(133) = me2(1,1);
   pars(134) = me2(1,2);
   pars(135) = me2(2,0);
   pars(136) = me2(2,1);
   pars(137) = me2(2,2);
   pars(138) = ms2;
   pars(139) = mH1I2(0,0);
   pars(140) = mH1I2(0,1);
   pars(141) = mH1I2(1,0);
   pars(142) = mH1I2(1,1);
   pars(143) = mH2I2(0,0);
   pars(144) = mH2I2(0,1);
   pars(145) = mH2I2(1,0);
   pars(146) = mH2I2(1,1);
   pars(147) = msI2(0,0);
   pars(148) = msI2(0,1);
   pars(149) = msI2(1,0);
   pars(150) = msI2(1,1);
   pars(151) = mDx2(0,0);
   pars(152) = mDx2(0,1);
   pars(153) = mDx2(0,2);
   pars(154) = mDx2(1,0);
   pars(155) = mDx2(1,1);
   pars(156) = mDx2(1,2);
   pars(157) = mDx2(2,0);
   pars(158) = mDx2(2,1);
   pars(159) = mDx2(2,2);
   pars(160) = mDxbar2(0,0);
   pars(161) = mDxbar2(0,1);
   pars(162) = mDxbar2(0,2);
   pars(163) = mDxbar2(1,0);
   pars(164) = mDxbar2(1,1);
   pars(165) = mDxbar2(1,2);
   pars(166) = mDxbar2(2,0);
   pars(167) = mDxbar2(2,1);
   pars(168) = mDxbar2(2,2);
   pars(169) = mHp2;
   pars(170) = mHpbar2;
   pars(171) = MassB;
   pars(172) = MassWB;
   pars(173) = MassG;
   pars(174) = MassBp;


   return pars;
}

void E6SSM_soft_parameters::print(std::ostream& ostr) const
{
   E6SSM_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "TYd = " << TYd << '\n';
   ostr << "TYe = " << TYe << '\n';
   ostr << "TKappa = " << TKappa << '\n';
   ostr << "TLambda12 = " << TLambda12 << '\n';
   ostr << "TLambdax = " << TLambdax << '\n';
   ostr << "TYu = " << TYu << '\n';
   ostr << "BMuPr = " << BMuPr << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "ms2 = " << ms2 << '\n';
   ostr << "mH1I2 = " << mH1I2 << '\n';
   ostr << "mH2I2 = " << mH2I2 << '\n';
   ostr << "msI2 = " << msI2 << '\n';
   ostr << "mDx2 = " << mDx2 << '\n';
   ostr << "mDxbar2 = " << mDxbar2 << '\n';
   ostr << "mHp2 = " << mHp2 << '\n';
   ostr << "mHpbar2 = " << mHpbar2 << '\n';
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "MassG = " << MassG << '\n';
   ostr << "MassBp = " << MassBp << '\n';

}

void E6SSM_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   E6SSM_susy_parameters::set(pars);

   TYd(0,0) = pars(49);
   TYd(0,1) = pars(50);
   TYd(0,2) = pars(51);
   TYd(1,0) = pars(52);
   TYd(1,1) = pars(53);
   TYd(1,2) = pars(54);
   TYd(2,0) = pars(55);
   TYd(2,1) = pars(56);
   TYd(2,2) = pars(57);
   TYe(0,0) = pars(58);
   TYe(0,1) = pars(59);
   TYe(0,2) = pars(60);
   TYe(1,0) = pars(61);
   TYe(1,1) = pars(62);
   TYe(1,2) = pars(63);
   TYe(2,0) = pars(64);
   TYe(2,1) = pars(65);
   TYe(2,2) = pars(66);
   TKappa(0,0) = pars(67);
   TKappa(0,1) = pars(68);
   TKappa(0,2) = pars(69);
   TKappa(1,0) = pars(70);
   TKappa(1,1) = pars(71);
   TKappa(1,2) = pars(72);
   TKappa(2,0) = pars(73);
   TKappa(2,1) = pars(74);
   TKappa(2,2) = pars(75);
   TLambda12(0,0) = pars(76);
   TLambda12(0,1) = pars(77);
   TLambda12(1,0) = pars(78);
   TLambda12(1,1) = pars(79);
   TLambdax = pars(80);
   TYu(0,0) = pars(81);
   TYu(0,1) = pars(82);
   TYu(0,2) = pars(83);
   TYu(1,0) = pars(84);
   TYu(1,1) = pars(85);
   TYu(1,2) = pars(86);
   TYu(2,0) = pars(87);
   TYu(2,1) = pars(88);
   TYu(2,2) = pars(89);
   BMuPr = pars(90);
   mq2(0,0) = pars(91);
   mq2(0,1) = pars(92);
   mq2(0,2) = pars(93);
   mq2(1,0) = pars(94);
   mq2(1,1) = pars(95);
   mq2(1,2) = pars(96);
   mq2(2,0) = pars(97);
   mq2(2,1) = pars(98);
   mq2(2,2) = pars(99);
   ml2(0,0) = pars(100);
   ml2(0,1) = pars(101);
   ml2(0,2) = pars(102);
   ml2(1,0) = pars(103);
   ml2(1,1) = pars(104);
   ml2(1,2) = pars(105);
   ml2(2,0) = pars(106);
   ml2(2,1) = pars(107);
   ml2(2,2) = pars(108);
   mHd2 = pars(109);
   mHu2 = pars(110);
   md2(0,0) = pars(111);
   md2(0,1) = pars(112);
   md2(0,2) = pars(113);
   md2(1,0) = pars(114);
   md2(1,1) = pars(115);
   md2(1,2) = pars(116);
   md2(2,0) = pars(117);
   md2(2,1) = pars(118);
   md2(2,2) = pars(119);
   mu2(0,0) = pars(120);
   mu2(0,1) = pars(121);
   mu2(0,2) = pars(122);
   mu2(1,0) = pars(123);
   mu2(1,1) = pars(124);
   mu2(1,2) = pars(125);
   mu2(2,0) = pars(126);
   mu2(2,1) = pars(127);
   mu2(2,2) = pars(128);
   me2(0,0) = pars(129);
   me2(0,1) = pars(130);
   me2(0,2) = pars(131);
   me2(1,0) = pars(132);
   me2(1,1) = pars(133);
   me2(1,2) = pars(134);
   me2(2,0) = pars(135);
   me2(2,1) = pars(136);
   me2(2,2) = pars(137);
   ms2 = pars(138);
   mH1I2(0,0) = pars(139);
   mH1I2(0,1) = pars(140);
   mH1I2(1,0) = pars(141);
   mH1I2(1,1) = pars(142);
   mH2I2(0,0) = pars(143);
   mH2I2(0,1) = pars(144);
   mH2I2(1,0) = pars(145);
   mH2I2(1,1) = pars(146);
   msI2(0,0) = pars(147);
   msI2(0,1) = pars(148);
   msI2(1,0) = pars(149);
   msI2(1,1) = pars(150);
   mDx2(0,0) = pars(151);
   mDx2(0,1) = pars(152);
   mDx2(0,2) = pars(153);
   mDx2(1,0) = pars(154);
   mDx2(1,1) = pars(155);
   mDx2(1,2) = pars(156);
   mDx2(2,0) = pars(157);
   mDx2(2,1) = pars(158);
   mDx2(2,2) = pars(159);
   mDxbar2(0,0) = pars(160);
   mDxbar2(0,1) = pars(161);
   mDxbar2(0,2) = pars(162);
   mDxbar2(1,0) = pars(163);
   mDxbar2(1,1) = pars(164);
   mDxbar2(1,2) = pars(165);
   mDxbar2(2,0) = pars(166);
   mDxbar2(2,1) = pars(167);
   mDxbar2(2,2) = pars(168);
   mHp2 = pars(169);
   mHpbar2 = pars(170);
   MassB = pars(171);
   MassWB = pars(172);
   MassG = pars(173);
   MassBp = pars(174);

}

E6SSM_soft_parameters::Soft_traces E6SSM_soft_parameters::calc_soft_traces(int loops) const
{
   Soft_traces soft_traces;

   if (loops > 0) {
      
      TRACE_STRUCT.Tr11 = Re(0.7745966692414834*g1*(-mHd2 - mHp2 + mHpbar2 + mHu2 + (
         md2).trace() - (mDx2).trace() + (mDxbar2).trace() + (me2).trace() - (mH1I2).
         trace() + (mH2I2).trace() - (ml2).trace() + (mq2).trace() - 2*(mu2).trace())
         );
      TRACE_STRUCT.Tr14 = Re(0.15811388300841897*gN*(-6*mHd2 + 4*mHp2 - 4*mHpbar2 - 4
         *mHu2 + 5*ms2 + 6*(md2).trace() - 6*(mDx2).trace() - 9*(mDxbar2).trace() + (
         me2).trace() - 6*(mH1I2).trace() - 4*(mH2I2).trace() + 4*(ml2).trace() + 6*(
         mq2).trace() + 5*(msI2).trace() + 3*(mu2).trace()));
      TRACE_STRUCT.Tr2U111 = Re(0.1*Sqr(g1)*(3*mHd2 + 3*mHp2 + 3*mHpbar2 + 3*mHu2 + 2
         *(md2).trace() + 2*(mDx2).trace() + 2*(mDxbar2).trace() + 6*(me2).trace() +
         3*(mH1I2).trace() + 3*(mH2I2).trace() + 3*(ml2).trace() + (mq2).trace() + 8*
         (mu2).trace()));
      TRACE_STRUCT.Tr2U114 = Re(0.1224744871391589*g1*gN*(3*mHd2 - 2*mHp2 - 2*mHpbar2
          - 2*mHu2 + 2*(md2).trace() + 2*(mDx2).trace() - 3*(mDxbar2).trace() + (me2)
         .trace() + 3*(mH1I2).trace() - 2*(mH2I2).trace() - 2*(ml2).trace() + (mq2).
         trace() - 2*(mu2).trace()));
      TRACE_STRUCT.Tr31 = Re(0.006454972243679028*g1*(60*(mHd2 - mHu2)*AbsSqr(Lambdax
         ) - 18*mHd2*Sqr(g1) - 18*mHp2*Sqr(g1) + 18*mHpbar2*Sqr(g1) + 18*mHu2*Sqr(g1)
         - 90*mHd2*Sqr(g2) - 90*mHp2*Sqr(g2) + 90*mHpbar2*Sqr(g2) + 90*mHu2*Sqr(g2) -
         27*mHd2*Sqr(gN) - 12*mHp2*Sqr(gN) + 12*mHpbar2*Sqr(gN) + 12*mHu2*Sqr(gN) + 4
         *(2*Sqr(g1) + 40*Sqr(g3) + 3*Sqr(gN))*(md2).trace() - 8*Sqr(g1)*(mDx2).trace
         () - 160*Sqr(g3)*(mDx2).trace() - 12*Sqr(gN)*(mDx2).trace() + 8*Sqr(g1)*(
         mDxbar2).trace() + 160*Sqr(g3)*(mDxbar2).trace() + 27*Sqr(gN)*(mDxbar2).
         trace() + 72*Sqr(g1)*(me2).trace() + 3*Sqr(gN)*(me2).trace() - 18*Sqr(g1)*(
         mH1I2).trace() - 90*Sqr(g2)*(mH1I2).trace() - 27*Sqr(gN)*(mH1I2).trace() +
         18*Sqr(g1)*(mH2I2).trace() + 90*Sqr(g2)*(mH2I2).trace() + 12*Sqr(gN)*(mH2I2)
         .trace() - 18*Sqr(g1)*(ml2).trace() - 90*Sqr(g2)*(ml2).trace() - 12*Sqr(gN)*
         (ml2).trace() + 2*Sqr(g1)*(mq2).trace() + 90*Sqr(g2)*(mq2).trace() + 160*Sqr
         (g3)*(mq2).trace() + 3*Sqr(gN)*(mq2).trace() - 64*Sqr(g1)*(mu2).trace() -
         320*Sqr(g3)*(mu2).trace() - 6*Sqr(gN)*(mu2).trace() + 180*mHd2*(Yd*Yd.
         adjoint()).trace() + 60*mHd2*(Ye*Ye.adjoint()).trace() - 180*mHu2*(Yu*Yu.
         adjoint()).trace() + 60*(mDx2*Kappa*(Kappa).adjoint()).trace() - 60*(mDxbar2
         *(Kappa).adjoint()*Kappa).trace() - 60*(mH2I2*Lambda12*(Lambda12).adjoint())
         .trace() - 120*(Yd*Yd.adjoint()*md2.conjugate()).trace() - 60*(Yd*mq2.
         conjugate()*Yd.adjoint()).trace() - 120*(Ye*Ye.adjoint()*me2.conjugate()).
         trace() + 60*(Ye*ml2.conjugate()*Ye.adjoint()).trace() + 240*(Yu*Yu.adjoint(
         )*mu2.conjugate()).trace() - 60*(Yu*mq2.conjugate()*Yu.adjoint()).trace() +
         60*(Lambda12*mH1I2.conjugate()*(Lambda12).adjoint()).trace()));
      TRACE_STRUCT.Tr22 = Re(0.5*(mHd2 + mHp2 + mHpbar2 + mHu2 + (mH1I2).trace() + (
         mH2I2).trace() + (ml2).trace() + 3*(mq2).trace()));
      TRACE_STRUCT.Tr23 = Re(0.5*((md2).trace() + (mDx2).trace() + (mDxbar2).trace()
         + 2*(mq2).trace() + (mu2).trace()));
      TRACE_STRUCT.Tr2U141 = Re(0.1224744871391589*g1*gN*(3*mHd2 - 2*mHp2 - 2*mHpbar2
          - 2*mHu2 + 2*(md2).trace() + 2*(mDx2).trace() - 3*(mDxbar2).trace() + (me2)
         .trace() + 3*(mH1I2).trace() - 2*(mH2I2).trace() - 2*(ml2).trace() + (mq2).
         trace() - 2*(mu2).trace()));
      TRACE_STRUCT.Tr2U144 = Re(0.025*Sqr(gN)*(18*mHd2 + 8*mHp2 + 8*mHpbar2 + 8*mHu2
         + 25*ms2 + 12*(md2).trace() + 12*(mDx2).trace() + 27*(mDxbar2).trace() + (
         me2).trace() + 18*(mH1I2).trace() + 8*(mH2I2).trace() + 8*(ml2).trace() + 6*
         (mq2).trace() + 25*(msI2).trace() + 3*(mu2).trace()));
      TRACE_STRUCT.Tr34 = Re(-0.003952847075210475*gN*(-40*(3*mHd2 + 2*mHu2 - 5*ms2)*
         AbsSqr(Lambdax) + 36*mHd2*Sqr(g1) - 24*mHp2*Sqr(g1) + 24*mHpbar2*Sqr(g1) +
         24*mHu2*Sqr(g1) + 180*mHd2*Sqr(g2) - 120*mHp2*Sqr(g2) + 120*mHpbar2*Sqr(g2)
         + 120*mHu2*Sqr(g2) + 54*mHd2*Sqr(gN) - 16*mHp2*Sqr(gN) + 16*mHpbar2*Sqr(gN)
         + 16*mHu2*Sqr(gN) - 125*ms2*Sqr(gN) - 8*(2*Sqr(g1) + 40*Sqr(g3) + 3*Sqr(gN))
         *(md2).trace() + 16*Sqr(g1)*(mDx2).trace() + 320*Sqr(g3)*(mDx2).trace() + 24
         *Sqr(gN)*(mDx2).trace() + 24*Sqr(g1)*(mDxbar2).trace() + 480*Sqr(g3)*(
         mDxbar2).trace() + 81*Sqr(gN)*(mDxbar2).trace() - 24*Sqr(g1)*(me2).trace() -
         Sqr(gN)*(me2).trace() + 36*Sqr(g1)*(mH1I2).trace() + 180*Sqr(g2)*(mH1I2).
         trace() + 54*Sqr(gN)*(mH1I2).trace() + 24*Sqr(g1)*(mH2I2).trace() + 120*Sqr(
         g2)*(mH2I2).trace() + 16*Sqr(gN)*(mH2I2).trace() - 24*Sqr(g1)*(ml2).trace()
         - 120*Sqr(g2)*(ml2).trace() - 16*Sqr(gN)*(ml2).trace() - 4*Sqr(g1)*(mq2).
         trace() - 180*Sqr(g2)*(mq2).trace() - 320*Sqr(g3)*(mq2).trace() - 6*Sqr(gN)*
         (mq2).trace() - 125*Sqr(gN)*(msI2).trace() - 32*Sqr(g1)*(mu2).trace() - 160*
         Sqr(g3)*(mu2).trace() - 3*Sqr(gN)*(mu2).trace() - 360*mHd2*(Yd*Yd.adjoint())
         .trace() - 120*mHd2*(Ye*Ye.adjoint()).trace() - 240*mHu2*(Yu*Yu.adjoint()).
         trace() + 300*ms2*(Kappa*(Kappa).adjoint()).trace() + 200*ms2*(Lambda12*(
         Lambda12).adjoint()).trace() - 120*(mDx2*Kappa*(Kappa).adjoint()).trace() -
         180*(mDxbar2*(Kappa).adjoint()*Kappa).trace() - 80*(mH2I2*Lambda12*(Lambda12
         ).adjoint()).trace() + 240*(Yd*Yd.adjoint()*md2.conjugate()).trace() + 120*(
         Yd*mq2.conjugate()*Yd.adjoint()).trace() + 40*(Ye*Ye.adjoint()*me2.conjugate
         ()).trace() + 80*(Ye*ml2.conjugate()*Ye.adjoint()).trace() + 120*(Yu*Yu.
         adjoint()*mu2.conjugate()).trace() + 120*(Yu*mq2.conjugate()*Yu.adjoint()).
         trace() - 120*(Lambda12*mH1I2.conjugate()*(Lambda12).adjoint()).trace()));

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceKappaAdjKappa = Re((Kappa*(Kappa).adjoint()).trace());
      TRACE_STRUCT.traceLambda12AdjLambda12 = Re((Lambda12*(Lambda12).adjoint()).
         trace());
      TRACE_STRUCT.traceAdjYdTYd = Re((Yd.adjoint()*TYd).trace());
      TRACE_STRUCT.traceAdjYeTYe = Re((Ye.adjoint()*TYe).trace());
      TRACE_STRUCT.traceAdjYuTYu = Re((Yu.adjoint()*TYu).trace());
      TRACE_STRUCT.traceAdjKappaTKappa = Re(((Kappa).adjoint()*TKappa).trace());
      TRACE_STRUCT.traceAdjLambda12TLambda12 = Re(((Lambda12).adjoint()*TLambda12).
         trace());
      TRACE_STRUCT.traceconjTYdTpTYd = Re((TYd.conjugate()*(TYd).transpose()).trace()
         );
      TRACE_STRUCT.traceconjTYeTpTYe = Re((TYe.conjugate()*(TYe).transpose()).trace()
         );
      TRACE_STRUCT.traceconjTYuTpTYu = Re((TYu.conjugate()*(TYu).transpose()).trace()
         );
      TRACE_STRUCT.traceconjTKappaTpTKappa = Re((TKappa.conjugate()*(TKappa).
         transpose()).trace());
      TRACE_STRUCT.traceconjTLambda12TpTLambda12 = Re((TLambda12.conjugate()*(
         TLambda12).transpose()).trace());
      TRACE_STRUCT.tracemd2YdAdjYd = Re((md2*Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceme2YeAdjYe = Re((me2*Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12 = Re((mH1I2*(Lambda12).adjoint()*
         Lambda12).trace());
      TRACE_STRUCT.traceml2AdjYeYe = Re((ml2*Ye.adjoint()*Ye).trace());
      TRACE_STRUCT.tracemq2AdjYdYd = Re((mq2*Yd.adjoint()*Yd).trace());
      TRACE_STRUCT.tracemq2AdjYuYu = Re((mq2*Yu.adjoint()*Yu).trace());
      TRACE_STRUCT.tracemu2YuAdjYu = Re((mu2*Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2 = Re((Kappa*(Kappa).adjoint()*mDx2.
         conjugate()).trace());
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa = Re((Kappa*mDxbar2.conjugate()*(
         Kappa).adjoint()).trace());
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2 = Re((Lambda12*(Lambda12).
         adjoint()*mH2I2.conjugate()).trace());

   }

   if (loops > 1) {
      TRACE_STRUCT.traceconjTYdTpYd = Re((TYd.conjugate()*Yd.transpose()).trace());
      TRACE_STRUCT.traceconjTYeTpYe = Re((TYe.conjugate()*Ye.transpose()).trace());
      TRACE_STRUCT.traceconjTYuTpYu = Re((TYu.conjugate()*Yu.transpose()).trace());
      TRACE_STRUCT.traceconjTKappaTpKappa = Re((TKappa.conjugate()*(Kappa).transpose(
         )).trace());
      TRACE_STRUCT.traceconjTLambda12TpLambda12 = Re((TLambda12.conjugate()*(Lambda12
         ).transpose()).trace());
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYdAdjYdTYdAdjYd = Re((Yd*Yd.adjoint()*TYd*Yd.adjoint()).trace
         ());
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd = Re((Yd*Yd.adjoint()*TYd*(TYd).adjoint()).
         trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYdAdjYuTYuAdjYd = Re((Yd*Yu.adjoint()*TYu*Yd.adjoint()).trace
         ());
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd = Re((Yd*Yu.adjoint()*TYu*(TYd).adjoint()).
         trace());
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd = Re((Yd*(TYd).adjoint()*TYd*Yd.adjoint()).
         trace());
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd = Re((Yd*(TYu).adjoint()*TYu*Yd.adjoint()).
         trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYeTYeAdjYe = Re((Ye*Ye.adjoint()*TYe*Ye.adjoint()).trace
         ());
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe = Re((Ye*Ye.adjoint()*TYe*(TYe).adjoint()).
         trace());
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe = Re((Ye*(TYe).adjoint()*TYe*Ye.adjoint()).
         trace());
      TRACE_STRUCT.traceYuAdjYdTYdAdjYu = Re((Yu*Yd.adjoint()*TYd*Yu.adjoint()).trace
         ());
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu = Re((Yu*Yd.adjoint()*TYd*(TYu).adjoint()).
         trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
         );
      TRACE_STRUCT.traceYuAdjYuTYuAdjYu = Re((Yu*Yu.adjoint()*TYu*Yu.adjoint()).trace
         ());
      TRACE_STRUCT.traceYuAdjYuTYuAdjTYu = Re((Yu*Yu.adjoint()*TYu*(TYu).adjoint()).
         trace());
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu = Re((Yu*(TYd).adjoint()*TYd*Yu.adjoint()).
         trace());
      TRACE_STRUCT.traceYuAdjTYuTYuAdjYu = Re((Yu*(TYu).adjoint()*TYu*Yu.adjoint()).
         trace());
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa = Re((Kappa*(Kappa).adjoint()*
         Kappa*(Kappa).adjoint()).trace());
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjKappa = Re((Kappa*(Kappa).adjoint()*
         TKappa*(Kappa).adjoint()).trace());
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjTKappa = Re((Kappa*(Kappa).adjoint()*
         TKappa*(TKappa).adjoint()).trace());
      TRACE_STRUCT.traceKappaAdjTKappaTKappaAdjKappa = Re((Kappa*(TKappa).adjoint()*
         TKappa*(Kappa).adjoint()).trace());
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12 = Re((Lambda12*(
         Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()).trace());
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjLambda12 = Re((Lambda12*(
         Lambda12).adjoint()*TLambda12*(Lambda12).adjoint()).trace());
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjTLambda12 = Re((Lambda12*(
         Lambda12).adjoint()*TLambda12*(TLambda12).adjoint()).trace());
      TRACE_STRUCT.traceLambda12AdjTLambda12TLambda12AdjLambda12 = Re((Lambda12*(
         TLambda12).adjoint()*TLambda12*(Lambda12).adjoint()).trace());
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd = Re((md2*Yd*Yd.adjoint()*Yd*Yd.adjoint()).
         trace());
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd = Re((md2*Yd*Yu.adjoint()*Yu*Yd.adjoint()).
         trace());
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe = Re((me2*Ye*Ye.adjoint()*Ye*Ye.adjoint()).
         trace());
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 = Re((mH1I2*(
         Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12).trace());
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe = Re((ml2*Ye.adjoint()*Ye*Ye.adjoint()*Ye).
         trace());
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd = Re((mq2*Yd.adjoint()*Yd*Yd.adjoint()*Yd).
         trace());
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu = Re((mq2*Yd.adjoint()*Yd*Yu.adjoint()*Yu).
         trace());
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd = Re((mq2*Yu.adjoint()*Yu*Yd.adjoint()*Yd).
         trace());
      TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu = Re((mq2*Yu.adjoint()*Yu*Yu.adjoint()*Yu).
         trace());
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu = Re((mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint()).
         trace());
      TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu = Re((mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint()).
         trace());
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappaconjmDx2 = Re((Kappa*(Kappa).
         adjoint()*Kappa*(Kappa).adjoint()*mDx2.conjugate()).trace());
      TRACE_STRUCT.traceKappaAdjKappaKappaconjmDxbar2AdjKappa = Re((Kappa*(Kappa).
         adjoint()*Kappa*mDxbar2.conjugate()*(Kappa).adjoint()).trace());
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2KappaAdjKappa = Re((Kappa*(Kappa).
         adjoint()*mDx2.conjugate()*Kappa*(Kappa).adjoint()).trace());
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappaKappaAdjKappa = Re((Kappa*mDxbar2.
         conjugate()*(Kappa).adjoint()*Kappa*(Kappa).adjoint()).trace());
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 = Re((
         Lambda12*(Lambda12).adjoint()*Lambda12*(Lambda12).adjoint()*mH2I2.conjugate(
         )).trace());
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 = Re((
         Lambda12*(Lambda12).adjoint()*mH2I2.conjugate()*Lambda12*(Lambda12).adjoint(
         )).trace());

   }

   if (loops > 2) {

   }

   return soft_traces;
}

std::ostream& operator<<(std::ostream& ostr, const E6SSM_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
