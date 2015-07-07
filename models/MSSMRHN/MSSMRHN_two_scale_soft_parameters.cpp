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

// File generated at Tue 7 Jul 2015 13:23:28

#include "MSSMRHN_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

MSSMRHN_soft_parameters::MSSMRHN_soft_parameters(const MSSMRHN_input_parameters& input_)
   : MSSMRHN_susy_parameters(input_)
   , TYd(Eigen::Matrix<double,3,3>::Zero()), TYe(Eigen::Matrix<double,3,3>
   ::Zero()), TYu(Eigen::Matrix<double,3,3>::Zero()), TYv(Eigen::Matrix<double,
   3,3>::Zero()), BMu(0), BMv(Eigen::Matrix<double,3,3>::Zero()), mq2(
   Eigen::Matrix<double,3,3>::Zero()), ml2(Eigen::Matrix<double,3,3>::Zero()),
   mHd2(0), mHu2(0), md2(Eigen::Matrix<double,3,3>::Zero()), mu2(Eigen::Matrix<
   double,3,3>::Zero()), me2(Eigen::Matrix<double,3,3>::Zero()), mv2(
   Eigen::Matrix<double,3,3>::Zero()), MassB(0), MassWB(0), MassG(0)

{
   set_number_of_parameters(numberOfParameters);
}

MSSMRHN_soft_parameters::MSSMRHN_soft_parameters(
   const MSSMRHN_susy_parameters& susy_model
   , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>&
   TYe_, const Eigen::Matrix<double,3,3>& TYu_, const Eigen::Matrix<double,3,3>
   & TYv_, double BMu_, const Eigen::Matrix<double,3,3>& BMv_, const
   Eigen::Matrix<double,3,3>& mq2_, const Eigen::Matrix<double,3,3>& ml2_,
   double mHd2_, double mHu2_, const Eigen::Matrix<double,3,3>& md2_, const
   Eigen::Matrix<double,3,3>& mu2_, const Eigen::Matrix<double,3,3>& me2_,
   const Eigen::Matrix<double,3,3>& mv2_, double MassB_, double MassWB_, double
   MassG_

)
   : MSSMRHN_susy_parameters(susy_model)
   , TYd(TYd_), TYe(TYe_), TYu(TYu_), TYv(TYv_), BMu(BMu_), BMv(BMv_), mq2(mq2_
   ), ml2(ml2_), mHd2(mHd2_), mHu2(mHu2_), md2(md2_), mu2(mu2_), me2(me2_), mv2
   (mv2_), MassB(MassB_), MassWB(MassWB_), MassG(MassG_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd MSSMRHN_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

MSSMRHN_soft_parameters MSSMRHN_soft_parameters::calc_beta() const
{
   Soft_traces soft_traces;
   calc_soft_traces(soft_traces);

   Eigen::Matrix<double,3,3> beta_TYd(calc_beta_TYd_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TYe(calc_beta_TYe_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TYu(calc_beta_TYu_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_TYv(calc_beta_TYv_one_loop(TRACE_STRUCT));
   double beta_BMu(calc_beta_BMu_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_BMv(calc_beta_BMv_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mq2(calc_beta_mq2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_ml2(calc_beta_ml2_one_loop(TRACE_STRUCT));
   double beta_mHd2(calc_beta_mHd2_one_loop(TRACE_STRUCT));
   double beta_mHu2(calc_beta_mHu2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_md2(calc_beta_md2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mu2(calc_beta_mu2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_me2(calc_beta_me2_one_loop(TRACE_STRUCT));
   Eigen::Matrix<double,3,3> beta_mv2(calc_beta_mv2_one_loop(TRACE_STRUCT));
   double beta_MassB(calc_beta_MassB_one_loop(TRACE_STRUCT));
   double beta_MassWB(calc_beta_MassWB_one_loop(TRACE_STRUCT));
   double beta_MassG(calc_beta_MassG_one_loop(TRACE_STRUCT));

   if (get_loops() > 1) {
      beta_TYd += calc_beta_TYd_two_loop(TRACE_STRUCT);
      beta_TYe += calc_beta_TYe_two_loop(TRACE_STRUCT);
      beta_TYu += calc_beta_TYu_two_loop(TRACE_STRUCT);
      beta_TYv += calc_beta_TYv_two_loop(TRACE_STRUCT);
      beta_BMu += calc_beta_BMu_two_loop(TRACE_STRUCT);
      beta_BMv += calc_beta_BMv_two_loop(TRACE_STRUCT);
      beta_mq2 += calc_beta_mq2_two_loop(TRACE_STRUCT);
      beta_ml2 += calc_beta_ml2_two_loop(TRACE_STRUCT);
      beta_mHd2 += calc_beta_mHd2_two_loop(TRACE_STRUCT);
      beta_mHu2 += calc_beta_mHu2_two_loop(TRACE_STRUCT);
      beta_md2 += calc_beta_md2_two_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_two_loop(TRACE_STRUCT);
      beta_me2 += calc_beta_me2_two_loop(TRACE_STRUCT);
      beta_mv2 += calc_beta_mv2_two_loop(TRACE_STRUCT);
      beta_MassB += calc_beta_MassB_two_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_two_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_two_loop(TRACE_STRUCT);

      if (get_loops() > 2) {

      }
   }


   const MSSMRHN_susy_parameters susy_betas(MSSMRHN_susy_parameters::calc_beta());

   return MSSMRHN_soft_parameters(susy_betas, beta_TYd, beta_TYe, beta_TYu, beta_TYv, beta_BMu, beta_BMv, beta_mq2, beta_ml2, beta_mHd2, beta_mHu2, beta_md2, beta_mu2, beta_me2, beta_mv2, beta_MassB, beta_MassWB, beta_MassG);
}

void MSSMRHN_soft_parameters::clear()
{
   MSSMRHN_susy_parameters::clear();

   TYd = Eigen::Matrix<double,3,3>::Zero();
   TYe = Eigen::Matrix<double,3,3>::Zero();
   TYu = Eigen::Matrix<double,3,3>::Zero();
   TYv = Eigen::Matrix<double,3,3>::Zero();
   BMu = 0.;
   BMv = Eigen::Matrix<double,3,3>::Zero();
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   mv2 = Eigen::Matrix<double,3,3>::Zero();
   MassB = 0.;
   MassWB = 0.;
   MassG = 0.;

}

const Eigen::ArrayXd MSSMRHN_soft_parameters::get() const
{
   Eigen::ArrayXd pars(MSSMRHN_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(51) = TYd(0,0);
   pars(52) = TYd(0,1);
   pars(53) = TYd(0,2);
   pars(54) = TYd(1,0);
   pars(55) = TYd(1,1);
   pars(56) = TYd(1,2);
   pars(57) = TYd(2,0);
   pars(58) = TYd(2,1);
   pars(59) = TYd(2,2);
   pars(60) = TYe(0,0);
   pars(61) = TYe(0,1);
   pars(62) = TYe(0,2);
   pars(63) = TYe(1,0);
   pars(64) = TYe(1,1);
   pars(65) = TYe(1,2);
   pars(66) = TYe(2,0);
   pars(67) = TYe(2,1);
   pars(68) = TYe(2,2);
   pars(69) = TYu(0,0);
   pars(70) = TYu(0,1);
   pars(71) = TYu(0,2);
   pars(72) = TYu(1,0);
   pars(73) = TYu(1,1);
   pars(74) = TYu(1,2);
   pars(75) = TYu(2,0);
   pars(76) = TYu(2,1);
   pars(77) = TYu(2,2);
   pars(78) = TYv(0,0);
   pars(79) = TYv(0,1);
   pars(80) = TYv(0,2);
   pars(81) = TYv(1,0);
   pars(82) = TYv(1,1);
   pars(83) = TYv(1,2);
   pars(84) = TYv(2,0);
   pars(85) = TYv(2,1);
   pars(86) = TYv(2,2);
   pars(87) = BMu;
   pars(88) = BMv(0,0);
   pars(89) = BMv(0,1);
   pars(90) = BMv(0,2);
   pars(91) = BMv(1,0);
   pars(92) = BMv(1,1);
   pars(93) = BMv(1,2);
   pars(94) = BMv(2,0);
   pars(95) = BMv(2,1);
   pars(96) = BMv(2,2);
   pars(97) = mq2(0,0);
   pars(98) = mq2(0,1);
   pars(99) = mq2(0,2);
   pars(100) = mq2(1,0);
   pars(101) = mq2(1,1);
   pars(102) = mq2(1,2);
   pars(103) = mq2(2,0);
   pars(104) = mq2(2,1);
   pars(105) = mq2(2,2);
   pars(106) = ml2(0,0);
   pars(107) = ml2(0,1);
   pars(108) = ml2(0,2);
   pars(109) = ml2(1,0);
   pars(110) = ml2(1,1);
   pars(111) = ml2(1,2);
   pars(112) = ml2(2,0);
   pars(113) = ml2(2,1);
   pars(114) = ml2(2,2);
   pars(115) = mHd2;
   pars(116) = mHu2;
   pars(117) = md2(0,0);
   pars(118) = md2(0,1);
   pars(119) = md2(0,2);
   pars(120) = md2(1,0);
   pars(121) = md2(1,1);
   pars(122) = md2(1,2);
   pars(123) = md2(2,0);
   pars(124) = md2(2,1);
   pars(125) = md2(2,2);
   pars(126) = mu2(0,0);
   pars(127) = mu2(0,1);
   pars(128) = mu2(0,2);
   pars(129) = mu2(1,0);
   pars(130) = mu2(1,1);
   pars(131) = mu2(1,2);
   pars(132) = mu2(2,0);
   pars(133) = mu2(2,1);
   pars(134) = mu2(2,2);
   pars(135) = me2(0,0);
   pars(136) = me2(0,1);
   pars(137) = me2(0,2);
   pars(138) = me2(1,0);
   pars(139) = me2(1,1);
   pars(140) = me2(1,2);
   pars(141) = me2(2,0);
   pars(142) = me2(2,1);
   pars(143) = me2(2,2);
   pars(144) = mv2(0,0);
   pars(145) = mv2(0,1);
   pars(146) = mv2(0,2);
   pars(147) = mv2(1,0);
   pars(148) = mv2(1,1);
   pars(149) = mv2(1,2);
   pars(150) = mv2(2,0);
   pars(151) = mv2(2,1);
   pars(152) = mv2(2,2);
   pars(153) = MassB;
   pars(154) = MassWB;
   pars(155) = MassG;


   return pars;
}

void MSSMRHN_soft_parameters::print(std::ostream& ostr) const
{
   MSSMRHN_susy_parameters::print(ostr);
   ostr << "soft parameters:\n";
   ostr << "TYd = " << TYd << '\n';
   ostr << "TYe = " << TYe << '\n';
   ostr << "TYu = " << TYu << '\n';
   ostr << "TYv = " << TYv << '\n';
   ostr << "BMu = " << BMu << '\n';
   ostr << "BMv = " << BMv << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "mv2 = " << mv2 << '\n';
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "MassG = " << MassG << '\n';

}

void MSSMRHN_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   MSSMRHN_susy_parameters::set(pars);

   TYd(0,0) = pars(51);
   TYd(0,1) = pars(52);
   TYd(0,2) = pars(53);
   TYd(1,0) = pars(54);
   TYd(1,1) = pars(55);
   TYd(1,2) = pars(56);
   TYd(2,0) = pars(57);
   TYd(2,1) = pars(58);
   TYd(2,2) = pars(59);
   TYe(0,0) = pars(60);
   TYe(0,1) = pars(61);
   TYe(0,2) = pars(62);
   TYe(1,0) = pars(63);
   TYe(1,1) = pars(64);
   TYe(1,2) = pars(65);
   TYe(2,0) = pars(66);
   TYe(2,1) = pars(67);
   TYe(2,2) = pars(68);
   TYu(0,0) = pars(69);
   TYu(0,1) = pars(70);
   TYu(0,2) = pars(71);
   TYu(1,0) = pars(72);
   TYu(1,1) = pars(73);
   TYu(1,2) = pars(74);
   TYu(2,0) = pars(75);
   TYu(2,1) = pars(76);
   TYu(2,2) = pars(77);
   TYv(0,0) = pars(78);
   TYv(0,1) = pars(79);
   TYv(0,2) = pars(80);
   TYv(1,0) = pars(81);
   TYv(1,1) = pars(82);
   TYv(1,2) = pars(83);
   TYv(2,0) = pars(84);
   TYv(2,1) = pars(85);
   TYv(2,2) = pars(86);
   BMu = pars(87);
   BMv(0,0) = pars(88);
   BMv(0,1) = pars(89);
   BMv(0,2) = pars(90);
   BMv(1,0) = pars(91);
   BMv(1,1) = pars(92);
   BMv(1,2) = pars(93);
   BMv(2,0) = pars(94);
   BMv(2,1) = pars(95);
   BMv(2,2) = pars(96);
   mq2(0,0) = pars(97);
   mq2(0,1) = pars(98);
   mq2(0,2) = pars(99);
   mq2(1,0) = pars(100);
   mq2(1,1) = pars(101);
   mq2(1,2) = pars(102);
   mq2(2,0) = pars(103);
   mq2(2,1) = pars(104);
   mq2(2,2) = pars(105);
   ml2(0,0) = pars(106);
   ml2(0,1) = pars(107);
   ml2(0,2) = pars(108);
   ml2(1,0) = pars(109);
   ml2(1,1) = pars(110);
   ml2(1,2) = pars(111);
   ml2(2,0) = pars(112);
   ml2(2,1) = pars(113);
   ml2(2,2) = pars(114);
   mHd2 = pars(115);
   mHu2 = pars(116);
   md2(0,0) = pars(117);
   md2(0,1) = pars(118);
   md2(0,2) = pars(119);
   md2(1,0) = pars(120);
   md2(1,1) = pars(121);
   md2(1,2) = pars(122);
   md2(2,0) = pars(123);
   md2(2,1) = pars(124);
   md2(2,2) = pars(125);
   mu2(0,0) = pars(126);
   mu2(0,1) = pars(127);
   mu2(0,2) = pars(128);
   mu2(1,0) = pars(129);
   mu2(1,1) = pars(130);
   mu2(1,2) = pars(131);
   mu2(2,0) = pars(132);
   mu2(2,1) = pars(133);
   mu2(2,2) = pars(134);
   me2(0,0) = pars(135);
   me2(0,1) = pars(136);
   me2(0,2) = pars(137);
   me2(1,0) = pars(138);
   me2(1,1) = pars(139);
   me2(1,2) = pars(140);
   me2(2,0) = pars(141);
   me2(2,1) = pars(142);
   me2(2,2) = pars(143);
   mv2(0,0) = pars(144);
   mv2(0,1) = pars(145);
   mv2(0,2) = pars(146);
   mv2(1,0) = pars(147);
   mv2(1,1) = pars(148);
   mv2(1,2) = pars(149);
   mv2(2,0) = pars(150);
   mv2(2,1) = pars(151);
   mv2(2,2) = pars(152);
   MassB = pars(153);
   MassWB = pars(154);
   MassG = pars(155);

}

void MSSMRHN_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
{
   TRACE_STRUCT.traceAdjYdTYd = Re((Yd.adjoint()*TYd).trace());
   TRACE_STRUCT.traceAdjYeTYe = Re((Ye.adjoint()*TYe).trace());
   TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
   TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
   TRACE_STRUCT.traceYdAdjYdTYdAdjYd = Re((Yd*Yd.adjoint()*TYd*Yd.adjoint())
      .trace());
   TRACE_STRUCT.traceYdAdjYuTYuAdjYd = Re((Yd*Yu.adjoint()*TYu*Yd.adjoint())
      .trace());
   TRACE_STRUCT.traceYeAdjYeTYeAdjYe = Re((Ye*Ye.adjoint()*TYe*Ye.adjoint())
      .trace());
   TRACE_STRUCT.traceYeAdjYvTYvAdjYe = Re((Ye*Yv.adjoint()*TYv*Ye.adjoint())
      .trace());
   TRACE_STRUCT.traceYuAdjYdTYdAdjYu = Re((Yu*Yd.adjoint()*TYd*Yu.adjoint())
      .trace());
   TRACE_STRUCT.traceYvAdjYeTYeAdjYv = Re((Yv*Ye.adjoint()*TYe*Yv.adjoint())
      .trace());
   TRACE_STRUCT.traceAdjYuTYu = Re((Yu.adjoint()*TYu).trace());
   TRACE_STRUCT.traceAdjYvTYv = Re((Yv.adjoint()*TYv).trace());
   TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
   TRACE_STRUCT.traceYvAdjYv = Re((Yv*Yv.adjoint()).trace());
   TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
      .trace());
   TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
      .trace());
   TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
      .trace());
   TRACE_STRUCT.traceYeAdjYvYvAdjYe = Re((Ye*Yv.adjoint()*Yv*Ye.adjoint())
      .trace());
   TRACE_STRUCT.traceYuAdjYuTYuAdjYu = Re((Yu*Yu.adjoint()*TYu*Yu.adjoint())
      .trace());
   TRACE_STRUCT.traceYvAdjYvTYvAdjYv = Re((Yv*Yv.adjoint()*TYv*Yv.adjoint())
      .trace());
   TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
      .trace());
   TRACE_STRUCT.traceYvAdjYvYvAdjYv = Re((Yv*Yv.adjoint()*Yv*Yv.adjoint())
      .trace());
   TRACE_STRUCT.traceconjTYdTpTYd = Re((TYd.conjugate()*(TYd).transpose())
      .trace());
   TRACE_STRUCT.traceconjTYeTpTYe = Re((TYe.conjugate()*(TYe).transpose())
      .trace());
   TRACE_STRUCT.tracemd2YdAdjYd = Re((md2*Yd*Yd.adjoint()).trace());
   TRACE_STRUCT.traceme2YeAdjYe = Re((me2*Ye*Ye.adjoint()).trace());
   TRACE_STRUCT.traceml2AdjYeYe = Re((ml2*Ye.adjoint()*Ye).trace());
   TRACE_STRUCT.tracemq2AdjYdYd = Re((mq2*Yd.adjoint()*Yd).trace());
   TRACE_STRUCT.traceconjTYdTpYd = Re((TYd.conjugate()*Yd.transpose()).trace())
      ;
   TRACE_STRUCT.traceconjTYeTpYe = Re((TYe.conjugate()*Ye.transpose()).trace())
      ;
   TRACE_STRUCT.traceconjTYuTpTYu = Re((TYu.conjugate()*(TYu).transpose())
      .trace());
   TRACE_STRUCT.traceconjTYvTpTYv = Re((TYv.conjugate()*(TYv).transpose())
      .trace());
   TRACE_STRUCT.traceml2AdjYvYv = Re((ml2*Yv.adjoint()*Yv).trace());
   TRACE_STRUCT.tracemq2AdjYuYu = Re((mq2*Yu.adjoint()*Yu).trace());
   TRACE_STRUCT.tracemu2YuAdjYu = Re((mu2*Yu*Yu.adjoint()).trace());
   TRACE_STRUCT.tracemv2YvAdjYv = Re((mv2*Yv*Yv.adjoint()).trace());
   TRACE_STRUCT.traceconjTYuTpYu = Re((TYu.conjugate()*Yu.transpose()).trace())
      ;
   TRACE_STRUCT.traceconjTYvTpYv = Re((TYv.conjugate()*Yv.transpose()).trace())
      ;
   TRACE_STRUCT.traceYdAdjYdTYdAdjTYd = Re((Yd*Yd.adjoint()*TYd*(TYd).adjoint()
      ).trace());
   TRACE_STRUCT.traceYdAdjYuTYuAdjTYd = Re((Yd*Yu.adjoint()*TYu*(TYd).adjoint()
      ).trace());
   TRACE_STRUCT.traceYdAdjTYdTYdAdjYd = Re((Yd*(TYd).adjoint()*TYd*Yd.adjoint()
      ).trace());
   TRACE_STRUCT.traceYdAdjTYuTYuAdjYd = Re((Yd*(TYu).adjoint()*TYu*Yd.adjoint()
      ).trace());
   TRACE_STRUCT.traceYeAdjYeTYeAdjTYe = Re((Ye*Ye.adjoint()*TYe*(TYe).adjoint()
      ).trace());
   TRACE_STRUCT.traceYeAdjYvTYvAdjTYe = Re((Ye*Yv.adjoint()*TYv*(TYe).adjoint()
      ).trace());
   TRACE_STRUCT.traceYeAdjTYeTYeAdjYe = Re((Ye*(TYe).adjoint()*TYe*Ye.adjoint()
      ).trace());
   TRACE_STRUCT.traceYeAdjTYvTYvAdjYe = Re((Ye*(TYv).adjoint()*TYv*Ye.adjoint()
      ).trace());
   TRACE_STRUCT.traceYuAdjYdTYdAdjTYu = Re((Yu*Yd.adjoint()*TYd*(TYu).adjoint()
      ).trace());
   TRACE_STRUCT.traceYuAdjTYdTYdAdjYu = Re((Yu*(TYd).adjoint()*TYd*Yu.adjoint()
      ).trace());
   TRACE_STRUCT.traceYvAdjYeTYeAdjTYv = Re((Yv*Ye.adjoint()*TYe*(TYv).adjoint()
      ).trace());
   TRACE_STRUCT.traceYvAdjTYeTYeAdjYv = Re((Yv*(TYe).adjoint()*TYe*Yv.adjoint()
      ).trace());
   TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd = Re((md2*Yd*Yd.adjoint()*Yd*Yd.adjoint(
      )).trace());
   TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd = Re((md2*Yd*Yu.adjoint()*Yu*Yd.adjoint(
      )).trace());
   TRACE_STRUCT.traceme2YeAdjYeYeAdjYe = Re((me2*Ye*Ye.adjoint()*Ye*Ye.adjoint(
      )).trace());
   TRACE_STRUCT.traceme2YeAdjYvYvAdjYe = Re((me2*Ye*Yv.adjoint()*Yv*Ye.adjoint(
      )).trace());
   TRACE_STRUCT.traceml2AdjYeYeAdjYeYe = Re((ml2*Ye.adjoint()*Ye*Ye.adjoint()*
      Ye).trace());
   TRACE_STRUCT.traceml2AdjYeYeAdjYvYv = Re((ml2*Ye.adjoint()*Ye*Yv.adjoint()*
      Yv).trace());
   TRACE_STRUCT.traceml2AdjYvYvAdjYeYe = Re((ml2*Yv.adjoint()*Yv*Ye.adjoint()*
      Ye).trace());
   TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd = Re((mq2*Yd.adjoint()*Yd*Yd.adjoint()*
      Yd).trace());
   TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu = Re((mq2*Yd.adjoint()*Yd*Yu.adjoint()*
      Yu).trace());
   TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd = Re((mq2*Yu.adjoint()*Yu*Yd.adjoint()*
      Yd).trace());
   TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu = Re((mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint(
      )).trace());
   TRACE_STRUCT.tracemv2YvAdjYeYeAdjYv = Re((mv2*Yv*Ye.adjoint()*Ye*Yv.adjoint(
      )).trace());
   TRACE_STRUCT.traceYuAdjYuTYuAdjTYu = Re((Yu*Yu.adjoint()*TYu*(TYu).adjoint()
      ).trace());
   TRACE_STRUCT.traceYuAdjTYuTYuAdjYu = Re((Yu*(TYu).adjoint()*TYu*Yu.adjoint()
      ).trace());
   TRACE_STRUCT.traceYvAdjYvTYvAdjTYv = Re((Yv*Yv.adjoint()*TYv*(TYv).adjoint()
      ).trace());
   TRACE_STRUCT.traceYvAdjTYvTYvAdjYv = Re((Yv*(TYv).adjoint()*TYv*Yv.adjoint()
      ).trace());
   TRACE_STRUCT.traceml2AdjYvYvAdjYvYv = Re((ml2*Yv.adjoint()*Yv*Yv.adjoint()*
      Yv).trace());
   TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu = Re((mq2*Yu.adjoint()*Yu*Yu.adjoint()*
      Yu).trace());
   TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu = Re((mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint(
      )).trace());
   TRACE_STRUCT.tracemv2YvAdjYvYvAdjYv = Re((mv2*Yv*Yv.adjoint()*Yv*Yv.adjoint(
      )).trace());


   TRACE_STRUCT.Tr11 = Re(0.7745966692414834*g1*(-mHd2 + mHu2 + (md2).trace() +
      (me2).trace() - (ml2).trace() + (mq2).trace() - 2*(mu2).trace()));
   TRACE_STRUCT.Tr2U111 = Re(0.1*Sqr(g1)*(3*mHd2 + 3*mHu2 + 2*(md2).trace() + 6
      *(me2).trace() + 3*(ml2).trace() + (mq2).trace() + 8*(mu2).trace()));
   TRACE_STRUCT.Tr31 = Re(0.012909944487358056*g1*(-9*mHd2*Sqr(g1) + 9*mHu2*Sqr
      (g1) - 45*mHd2*Sqr(g2) + 45*mHu2*Sqr(g2) + 4*(Sqr(g1) + 20*Sqr(g3))*(md2)
      .trace() + 36*Sqr(g1)*(me2).trace() - 9*Sqr(g1)*(ml2).trace() - 45*Sqr(g2)*(
      ml2).trace() + Sqr(g1)*(mq2).trace() + 45*Sqr(g2)*(mq2).trace() + 80*Sqr(g3)
      *(mq2).trace() - 32*Sqr(g1)*(mu2).trace() - 160*Sqr(g3)*(mu2).trace() + 90*
      mHd2*(Yd*Yd.adjoint()).trace() + 30*mHd2*(Ye*Ye.adjoint()).trace() - 90*mHu2
      *(Yu*Yu.adjoint()).trace() - 30*mHu2*(Yv*Yv.adjoint()).trace() - 60*(Yd*
      Yd.adjoint()*md2.conjugate()).trace() - 30*(Yd*mq2.conjugate()*Yd.adjoint())
      .trace() - 60*(Ye*Ye.adjoint()*me2.conjugate()).trace() + 30*(Ye*
      ml2.conjugate()*Ye.adjoint()).trace() + 120*(Yu*Yu.adjoint()*mu2.conjugate()
      ).trace() - 30*(Yu*mq2.conjugate()*Yu.adjoint()).trace() + 30*(Yv*
      ml2.conjugate()*Yv.adjoint()).trace()));
   TRACE_STRUCT.Tr22 = Re(0.5*(mHd2 + mHu2 + (ml2).trace() + 3*(mq2).trace()));
   TRACE_STRUCT.Tr23 = Re(0.5*((md2).trace() + 2*(mq2).trace() + (mu2).trace())
      );

}

std::ostream& operator<<(std::ostream& ostr, const MSSMRHN_soft_parameters& soft_pars)
{
   soft_pars.print(std::cout);
   return ostr;
}

} // namespace flexiblesusy
