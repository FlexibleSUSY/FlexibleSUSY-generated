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

// File generated at Mon 9 May 2016 12:42:32

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES() calc_soft_traces(TRACE_STRUCT);

const int UMSSM_soft_parameters::numberOfParameters;

UMSSM_soft_parameters::UMSSM_soft_parameters(const UMSSM_input_parameters& input_)
   : UMSSM_susy_parameters(input_)
   , TYd(Eigen::Matrix<double,3,3>::Zero()), TYe(Eigen::Matrix<double,3,3>
   ::Zero()), TLambdax(0), TYv(Eigen::Matrix<double,3,3>::Zero()), TYu(
   Eigen::Matrix<double,3,3>::Zero()), mq2(Eigen::Matrix<double,3,3>::Zero()),
   ml2(Eigen::Matrix<double,3,3>::Zero()), mHd2(0), mHu2(0), md2(Eigen::Matrix<
   double,3,3>::Zero()), mu2(Eigen::Matrix<double,3,3>::Zero()), me2(
   Eigen::Matrix<double,3,3>::Zero()), mvR2(Eigen::Matrix<double,3,3>::Zero()),
   ms2(0), MassB(0), MassWB(0), MassG(0), MassU(0)

{
   set_number_of_parameters(numberOfParameters);
}

UMSSM_soft_parameters::UMSSM_soft_parameters(
   const UMSSM_susy_parameters& susy_model
   , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>&
   TYe_, double TLambdax_, const Eigen::Matrix<double,3,3>& TYv_, const
   Eigen::Matrix<double,3,3>& TYu_, const Eigen::Matrix<double,3,3>& mq2_,
   const Eigen::Matrix<double,3,3>& ml2_, double mHd2_, double mHu2_, const
   Eigen::Matrix<double,3,3>& md2_, const Eigen::Matrix<double,3,3>& mu2_,
   const Eigen::Matrix<double,3,3>& me2_, const Eigen::Matrix<double,3,3>&
   mvR2_, double ms2_, double MassB_, double MassWB_, double MassG_, double
   MassU_

)
   : UMSSM_susy_parameters(susy_model)
   , TYd(TYd_), TYe(TYe_), TLambdax(TLambdax_), TYv(TYv_), TYu(TYu_), mq2(mq2_)
   , ml2(ml2_), mHd2(mHd2_), mHu2(mHu2_), md2(md2_), mu2(mu2_), me2(me2_), mvR2
   (mvR2_), ms2(ms2_), MassB(MassB_), MassWB(MassWB_), MassG(MassG_), MassU(
   MassU_)

{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd UMSSM_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

UMSSM_soft_parameters UMSSM_soft_parameters::calc_beta() const
{
   Eigen::Matrix<double,3,3> beta_TYd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_TYe = Eigen::Matrix<double,3,3>::Zero();
   double beta_TLambdax = 0.;
   Eigen::Matrix<double,3,3> beta_TYv = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_TYu = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_mq2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_ml2 = Eigen::Matrix<double,3,3>::Zero();
   double beta_mHd2 = 0.;
   double beta_mHu2 = 0.;
   Eigen::Matrix<double,3,3> beta_md2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_mu2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_me2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_mvR2 = Eigen::Matrix<double,3,3>::Zero();
   double beta_ms2 = 0.;
   double beta_MassB = 0.;
   double beta_MassWB = 0.;
   double beta_MassG = 0.;
   double beta_MassU = 0.;

   if (get_loops() > 0) {
      TRACE_STRUCT_TYPE TRACE_STRUCT;
      CALCULATE_TRACES();

      beta_TYd += calc_beta_TYd_one_loop(TRACE_STRUCT);
      beta_TYe += calc_beta_TYe_one_loop(TRACE_STRUCT);
      beta_TLambdax += calc_beta_TLambdax_one_loop(TRACE_STRUCT);
      beta_TYv += calc_beta_TYv_one_loop(TRACE_STRUCT);
      beta_TYu += calc_beta_TYu_one_loop(TRACE_STRUCT);
      beta_mq2 += calc_beta_mq2_one_loop(TRACE_STRUCT);
      beta_ml2 += calc_beta_ml2_one_loop(TRACE_STRUCT);
      beta_mHd2 += calc_beta_mHd2_one_loop(TRACE_STRUCT);
      beta_mHu2 += calc_beta_mHu2_one_loop(TRACE_STRUCT);
      beta_md2 += calc_beta_md2_one_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_one_loop(TRACE_STRUCT);
      beta_me2 += calc_beta_me2_one_loop(TRACE_STRUCT);
      beta_mvR2 += calc_beta_mvR2_one_loop(TRACE_STRUCT);
      beta_ms2 += calc_beta_ms2_one_loop(TRACE_STRUCT);
      beta_MassB += calc_beta_MassB_one_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_one_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_one_loop(TRACE_STRUCT);
      beta_MassU += calc_beta_MassU_one_loop(TRACE_STRUCT);

      if (get_loops() > 1) {
         beta_TYd += calc_beta_TYd_two_loop(TRACE_STRUCT);
         beta_TYe += calc_beta_TYe_two_loop(TRACE_STRUCT);
         beta_TLambdax += calc_beta_TLambdax_two_loop(TRACE_STRUCT);
         beta_TYv += calc_beta_TYv_two_loop(TRACE_STRUCT);
         beta_TYu += calc_beta_TYu_two_loop(TRACE_STRUCT);
         beta_mq2 += calc_beta_mq2_two_loop(TRACE_STRUCT);
         beta_ml2 += calc_beta_ml2_two_loop(TRACE_STRUCT);
         beta_mHd2 += calc_beta_mHd2_two_loop(TRACE_STRUCT);
         beta_mHu2 += calc_beta_mHu2_two_loop(TRACE_STRUCT);
         beta_md2 += calc_beta_md2_two_loop(TRACE_STRUCT);
         beta_mu2 += calc_beta_mu2_two_loop(TRACE_STRUCT);
         beta_me2 += calc_beta_me2_two_loop(TRACE_STRUCT);
         beta_mvR2 += calc_beta_mvR2_two_loop(TRACE_STRUCT);
         beta_ms2 += calc_beta_ms2_two_loop(TRACE_STRUCT);
         beta_MassB += calc_beta_MassB_two_loop(TRACE_STRUCT);
         beta_MassWB += calc_beta_MassWB_two_loop(TRACE_STRUCT);
         beta_MassG += calc_beta_MassG_two_loop(TRACE_STRUCT);
         beta_MassU += calc_beta_MassU_two_loop(TRACE_STRUCT);

         if (get_loops() > 2) {

         }
      }
   }


   const UMSSM_susy_parameters susy_betas(UMSSM_susy_parameters::calc_beta());

   return UMSSM_soft_parameters(susy_betas, beta_TYd, beta_TYe, beta_TLambdax, beta_TYv, beta_TYu, beta_mq2, beta_ml2, beta_mHd2, beta_mHu2, beta_md2, beta_mu2, beta_me2, beta_mvR2, beta_ms2, beta_MassB, beta_MassWB, beta_MassG, beta_MassU);
}

void UMSSM_soft_parameters::clear()
{
   UMSSM_susy_parameters::clear();

   TYd = Eigen::Matrix<double,3,3>::Zero();
   TYe = Eigen::Matrix<double,3,3>::Zero();
   TLambdax = 0.;
   TYv = Eigen::Matrix<double,3,3>::Zero();
   TYu = Eigen::Matrix<double,3,3>::Zero();
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   mvR2 = Eigen::Matrix<double,3,3>::Zero();
   ms2 = 0.;
   MassB = 0.;
   MassWB = 0.;
   MassG = 0.;
   MassU = 0.;

}

Eigen::ArrayXd UMSSM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(UMSSM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(44) = TYd(0,0);
   pars(45) = TYd(0,1);
   pars(46) = TYd(0,2);
   pars(47) = TYd(1,0);
   pars(48) = TYd(1,1);
   pars(49) = TYd(1,2);
   pars(50) = TYd(2,0);
   pars(51) = TYd(2,1);
   pars(52) = TYd(2,2);
   pars(53) = TYe(0,0);
   pars(54) = TYe(0,1);
   pars(55) = TYe(0,2);
   pars(56) = TYe(1,0);
   pars(57) = TYe(1,1);
   pars(58) = TYe(1,2);
   pars(59) = TYe(2,0);
   pars(60) = TYe(2,1);
   pars(61) = TYe(2,2);
   pars(62) = TLambdax;
   pars(63) = TYv(0,0);
   pars(64) = TYv(0,1);
   pars(65) = TYv(0,2);
   pars(66) = TYv(1,0);
   pars(67) = TYv(1,1);
   pars(68) = TYv(1,2);
   pars(69) = TYv(2,0);
   pars(70) = TYv(2,1);
   pars(71) = TYv(2,2);
   pars(72) = TYu(0,0);
   pars(73) = TYu(0,1);
   pars(74) = TYu(0,2);
   pars(75) = TYu(1,0);
   pars(76) = TYu(1,1);
   pars(77) = TYu(1,2);
   pars(78) = TYu(2,0);
   pars(79) = TYu(2,1);
   pars(80) = TYu(2,2);
   pars(81) = mq2(0,0);
   pars(82) = mq2(0,1);
   pars(83) = mq2(0,2);
   pars(84) = mq2(1,0);
   pars(85) = mq2(1,1);
   pars(86) = mq2(1,2);
   pars(87) = mq2(2,0);
   pars(88) = mq2(2,1);
   pars(89) = mq2(2,2);
   pars(90) = ml2(0,0);
   pars(91) = ml2(0,1);
   pars(92) = ml2(0,2);
   pars(93) = ml2(1,0);
   pars(94) = ml2(1,1);
   pars(95) = ml2(1,2);
   pars(96) = ml2(2,0);
   pars(97) = ml2(2,1);
   pars(98) = ml2(2,2);
   pars(99) = mHd2;
   pars(100) = mHu2;
   pars(101) = md2(0,0);
   pars(102) = md2(0,1);
   pars(103) = md2(0,2);
   pars(104) = md2(1,0);
   pars(105) = md2(1,1);
   pars(106) = md2(1,2);
   pars(107) = md2(2,0);
   pars(108) = md2(2,1);
   pars(109) = md2(2,2);
   pars(110) = mu2(0,0);
   pars(111) = mu2(0,1);
   pars(112) = mu2(0,2);
   pars(113) = mu2(1,0);
   pars(114) = mu2(1,1);
   pars(115) = mu2(1,2);
   pars(116) = mu2(2,0);
   pars(117) = mu2(2,1);
   pars(118) = mu2(2,2);
   pars(119) = me2(0,0);
   pars(120) = me2(0,1);
   pars(121) = me2(0,2);
   pars(122) = me2(1,0);
   pars(123) = me2(1,1);
   pars(124) = me2(1,2);
   pars(125) = me2(2,0);
   pars(126) = me2(2,1);
   pars(127) = me2(2,2);
   pars(128) = mvR2(0,0);
   pars(129) = mvR2(0,1);
   pars(130) = mvR2(0,2);
   pars(131) = mvR2(1,0);
   pars(132) = mvR2(1,1);
   pars(133) = mvR2(1,2);
   pars(134) = mvR2(2,0);
   pars(135) = mvR2(2,1);
   pars(136) = mvR2(2,2);
   pars(137) = ms2;
   pars(138) = MassB;
   pars(139) = MassWB;
   pars(140) = MassG;
   pars(141) = MassU;


   return pars;
}

void UMSSM_soft_parameters::print(std::ostream& ostr) const
{
   UMSSM_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "TYd = " << TYd << '\n';
   ostr << "TYe = " << TYe << '\n';
   ostr << "TLambdax = " << TLambdax << '\n';
   ostr << "TYv = " << TYv << '\n';
   ostr << "TYu = " << TYu << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "mvR2 = " << mvR2 << '\n';
   ostr << "ms2 = " << ms2 << '\n';
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "MassG = " << MassG << '\n';
   ostr << "MassU = " << MassU << '\n';

}

void UMSSM_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   UMSSM_susy_parameters::set(pars);

   TYd(0,0) = pars(44);
   TYd(0,1) = pars(45);
   TYd(0,2) = pars(46);
   TYd(1,0) = pars(47);
   TYd(1,1) = pars(48);
   TYd(1,2) = pars(49);
   TYd(2,0) = pars(50);
   TYd(2,1) = pars(51);
   TYd(2,2) = pars(52);
   TYe(0,0) = pars(53);
   TYe(0,1) = pars(54);
   TYe(0,2) = pars(55);
   TYe(1,0) = pars(56);
   TYe(1,1) = pars(57);
   TYe(1,2) = pars(58);
   TYe(2,0) = pars(59);
   TYe(2,1) = pars(60);
   TYe(2,2) = pars(61);
   TLambdax = pars(62);
   TYv(0,0) = pars(63);
   TYv(0,1) = pars(64);
   TYv(0,2) = pars(65);
   TYv(1,0) = pars(66);
   TYv(1,1) = pars(67);
   TYv(1,2) = pars(68);
   TYv(2,0) = pars(69);
   TYv(2,1) = pars(70);
   TYv(2,2) = pars(71);
   TYu(0,0) = pars(72);
   TYu(0,1) = pars(73);
   TYu(0,2) = pars(74);
   TYu(1,0) = pars(75);
   TYu(1,1) = pars(76);
   TYu(1,2) = pars(77);
   TYu(2,0) = pars(78);
   TYu(2,1) = pars(79);
   TYu(2,2) = pars(80);
   mq2(0,0) = pars(81);
   mq2(0,1) = pars(82);
   mq2(0,2) = pars(83);
   mq2(1,0) = pars(84);
   mq2(1,1) = pars(85);
   mq2(1,2) = pars(86);
   mq2(2,0) = pars(87);
   mq2(2,1) = pars(88);
   mq2(2,2) = pars(89);
   ml2(0,0) = pars(90);
   ml2(0,1) = pars(91);
   ml2(0,2) = pars(92);
   ml2(1,0) = pars(93);
   ml2(1,1) = pars(94);
   ml2(1,2) = pars(95);
   ml2(2,0) = pars(96);
   ml2(2,1) = pars(97);
   ml2(2,2) = pars(98);
   mHd2 = pars(99);
   mHu2 = pars(100);
   md2(0,0) = pars(101);
   md2(0,1) = pars(102);
   md2(0,2) = pars(103);
   md2(1,0) = pars(104);
   md2(1,1) = pars(105);
   md2(1,2) = pars(106);
   md2(2,0) = pars(107);
   md2(2,1) = pars(108);
   md2(2,2) = pars(109);
   mu2(0,0) = pars(110);
   mu2(0,1) = pars(111);
   mu2(0,2) = pars(112);
   mu2(1,0) = pars(113);
   mu2(1,1) = pars(114);
   mu2(1,2) = pars(115);
   mu2(2,0) = pars(116);
   mu2(2,1) = pars(117);
   mu2(2,2) = pars(118);
   me2(0,0) = pars(119);
   me2(0,1) = pars(120);
   me2(0,2) = pars(121);
   me2(1,0) = pars(122);
   me2(1,1) = pars(123);
   me2(1,2) = pars(124);
   me2(2,0) = pars(125);
   me2(2,1) = pars(126);
   me2(2,2) = pars(127);
   mvR2(0,0) = pars(128);
   mvR2(0,1) = pars(129);
   mvR2(0,2) = pars(130);
   mvR2(1,0) = pars(131);
   mvR2(1,1) = pars(132);
   mvR2(1,2) = pars(133);
   mvR2(2,0) = pars(134);
   mvR2(2,1) = pars(135);
   mvR2(2,2) = pars(136);
   ms2 = pars(137);
   MassB = pars(138);
   MassWB = pars(139);
   MassG = pars(140);
   MassU = pars(141);

}

void UMSSM_soft_parameters::calc_soft_traces(Soft_traces& soft_traces) const
{
   if (get_loops() > 0) {
      const auto QHd = INPUT(QHd);
      const auto QHu = INPUT(QHu);
      const auto Qs = INPUT(Qs);
      const auto Qd = INPUT(Qd);
      const auto Qe = INPUT(Qe);
      const auto Ql = INPUT(Ql);
      const auto Qq = INPUT(Qq);
      const auto Qu = INPUT(Qu);
      const auto Qv = INPUT(Qv);

      TRACE_STRUCT.Tr11 = Re(0.7745966692414834*g1*(-mHd2 + mHu2 + (md2).trace() +
         (me2).trace() - (ml2).trace() + (mq2).trace() - 2*(mu2).trace()));
      TRACE_STRUCT.Tr14 = Re(gp*(2*mHd2*QHd + 2*mHu2*QHu + ms2*Qs + 3*Qd*(md2)
         .trace() + Qe*(me2).trace() + 2*Ql*(ml2).trace() + 6*Qq*(mq2).trace() + 3*Qu
         *(mu2).trace() + Qv*(mvR2).trace()));
      TRACE_STRUCT.Tr2U111 = Re(0.1*Sqr(g1)*(3*mHd2 + 3*mHu2 + 2*(md2).trace() + 6
         *(me2).trace() + 3*(ml2).trace() + (mq2).trace() + 8*(mu2).trace()));
      TRACE_STRUCT.Tr2U114 = Re(0.7745966692414834*g1*gp*(-(mHd2*QHd) + mHu2*QHu +
         Qd*(md2).trace() + Qe*(me2).trace() - Ql*(ml2).trace() + Qq*(mq2).trace() -
         2*Qu*(mu2).trace()));
      TRACE_STRUCT.Tr31 = Re(0.012909944487358056*g1*(30*(mHd2 - mHu2)*AbsSqr(
         Lambdax) - 9*mHd2*Sqr(g1) + 9*mHu2*Sqr(g1) - 45*mHd2*Sqr(g2) + 45*mHu2*Sqr(
         g2) - 60*mHd2*Sqr(gp)*Sqr(QHd) + 60*mHu2*Sqr(gp)*Sqr(QHu) + 4*(Sqr(g1) + 20*
         Sqr(g3) + 15*Sqr(gp)*Sqr(Qd))*(md2).trace() + 36*Sqr(g1)*(me2).trace() + 60*
         Sqr(gp)*Sqr(Qe)*(me2).trace() - 9*Sqr(g1)*(ml2).trace() - 45*Sqr(g2)*(ml2)
         .trace() - 60*Sqr(gp)*Sqr(Ql)*(ml2).trace() + Sqr(g1)*(mq2).trace() + 45*Sqr
         (g2)*(mq2).trace() + 80*Sqr(g3)*(mq2).trace() + 60*Sqr(gp)*Sqr(Qq)*(mq2)
         .trace() - 32*Sqr(g1)*(mu2).trace() - 160*Sqr(g3)*(mu2).trace() - 120*Sqr(gp
         )*Sqr(Qu)*(mu2).trace() + 90*mHd2*(Yd*Yd.adjoint()).trace() + 30*mHd2*(Ye*
         Ye.adjoint()).trace() - 90*mHu2*(Yu*Yu.adjoint()).trace() - 30*mHu2*(Yv*
         Yv.adjoint()).trace() + 30*(ml2*Yv*Yv.adjoint()).trace() - 60*(Yd*Yd.adjoint
         ()*md2.conjugate()).trace() - 30*(Yd*mq2.conjugate()*Yd.adjoint()).trace() -
         60*(Ye*Ye.adjoint()*me2.conjugate()).trace() + 30*(Ye*ml2.conjugate()*
         Ye.adjoint()).trace() + 120*(Yu*Yu.adjoint()*mu2.conjugate()).trace() - 30*(
         Yu*mq2.conjugate()*Yu.adjoint()).trace()));
      TRACE_STRUCT.Tr22 = Re(0.5*(mHd2 + mHu2 + (ml2).trace() + 3*(mq2).trace()));
      TRACE_STRUCT.Tr23 = Re(0.5*((md2).trace() + 2*(mq2).trace() + (mu2).trace())
         );
      TRACE_STRUCT.Tr2U141 = Re(0.7745966692414834*g1*gp*(-(mHd2*QHd) + mHu2*QHu +
         Qd*(md2).trace() + Qe*(me2).trace() - Ql*(ml2).trace() + Qq*(mq2).trace() -
         2*Qu*(mu2).trace()));
      TRACE_STRUCT.Tr2U144 = Re(Sqr(gp)*(2*mHd2*Sqr(QHd) + 2*mHu2*Sqr(QHu) + ms2*
         Sqr(Qs) + 3*Sqr(Qd)*(md2).trace() + Sqr(Qe)*(me2).trace() + 2*Sqr(Ql)*(ml2)
         .trace() + 6*Sqr(Qq)*(mq2).trace() + 3*Sqr(Qu)*(mu2).trace() + Sqr(Qv)*(mvR2
         ).trace()));
      TRACE_STRUCT.Tr34 = Re(0.1*gp*(-10*(mHd2*QHd + mHu2*QHu + ms2*Qs)*AbsSqr(
         Lambdax) + 3*mHd2*QHd*Sqr(g1) + 3*mHu2*QHu*Sqr(g1) + 15*mHd2*QHd*Sqr(g2) +
         15*mHu2*QHu*Sqr(g2) + 20*mHd2*Power(QHd,3)*Sqr(gp) + 20*mHu2*Power(QHu,3)*
         Sqr(gp) + 10*ms2*Power(Qs,3)*Sqr(gp) + 2*Qd*(Sqr(g1) + 20*Sqr(g3) + 15*Sqr(
         gp)*Sqr(Qd))*(md2).trace() + 6*Qe*Sqr(g1)*(me2).trace() + 10*Power(Qe,3)*Sqr
         (gp)*(me2).trace() + 3*Ql*Sqr(g1)*(ml2).trace() + 15*Ql*Sqr(g2)*(ml2).trace(
         ) + 20*Power(Ql,3)*Sqr(gp)*(ml2).trace() + Qq*Sqr(g1)*(mq2).trace() + 45*Qq*
         Sqr(g2)*(mq2).trace() + 80*Qq*Sqr(g3)*(mq2).trace() + 60*Power(Qq,3)*Sqr(gp)
         *(mq2).trace() + 8*Qu*Sqr(g1)*(mu2).trace() + 40*Qu*Sqr(g3)*(mu2).trace() +
         30*Power(Qu,3)*Sqr(gp)*(mu2).trace() + 10*Power(Qv,3)*Sqr(gp)*(mvR2).trace()
         - 30*mHd2*QHd*(Yd*Yd.adjoint()).trace() - 10*mHd2*QHd*(Ye*Ye.adjoint())
         .trace() - 30*mHu2*QHu*(Yu*Yu.adjoint()).trace() - 10*mHu2*QHu*(Yv*
         Yv.adjoint()).trace() - 10*Ql*(ml2*Yv*Yv.adjoint()).trace() - 10*Qv*(mvR2*
         Yv.adjoint()*Yv).trace() - 30*Qd*(Yd*Yd.adjoint()*md2.conjugate()).trace() -
         30*Qq*(Yd*mq2.conjugate()*Yd.adjoint()).trace() - 10*Qe*(Ye*Ye.adjoint()*
         me2.conjugate()).trace() - 10*Ql*(Ye*ml2.conjugate()*Ye.adjoint()).trace() -
         30*Qu*(Yu*Yu.adjoint()*mu2.conjugate()).trace() - 30*Qq*(Yu*mq2.conjugate()
         *Yu.adjoint()).trace()));

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceYvAdjYv = Re((Yv*Yv.adjoint()).trace());
      TRACE_STRUCT.traceAdjYdTYd = Re((Yd.adjoint()*TYd).trace());
      TRACE_STRUCT.traceAdjYeTYe = Re((Ye.adjoint()*TYe).trace());
      TRACE_STRUCT.traceAdjYuTYu = Re((Yu.adjoint()*TYu).trace());
      TRACE_STRUCT.traceAdjYvTYv = Re((Yv.adjoint()*TYv).trace());
      TRACE_STRUCT.traceconjTYdTpTYd = Re((TYd.conjugate()*(TYd).transpose())
         .trace());
      TRACE_STRUCT.traceconjTYeTpTYe = Re((TYe.conjugate()*(TYe).transpose())
         .trace());
      TRACE_STRUCT.traceconjTYuTpTYu = Re((TYu.conjugate()*(TYu).transpose())
         .trace());
      TRACE_STRUCT.traceconjTYvTpTYv = Re((TYv.conjugate()*(TYv).transpose())
         .trace());
      TRACE_STRUCT.tracemd2YdAdjYd = Re((md2*Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceme2YeAdjYe = Re((me2*Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceml2AdjYeYe = Re((ml2*Ye.adjoint()*Ye).trace());
      TRACE_STRUCT.tracemq2AdjYdYd = Re((mq2*Yd.adjoint()*Yd).trace());
      TRACE_STRUCT.tracemq2AdjYuYu = Re((mq2*Yu.adjoint()*Yu).trace());
      TRACE_STRUCT.tracemu2YuAdjYu = Re((mu2*Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceYvAdjYvconjml2 = Re((Yv*Yv.adjoint()*ml2.conjugate())
         .trace());
      TRACE_STRUCT.traceYvconjmvR2AdjYv = Re((Yv*mvR2.conjugate()*Yv.adjoint())
         .trace());

   }

   if (get_loops() > 1) {
      TRACE_STRUCT.traceconjTYdTpYd = Re((TYd.conjugate()*Yd.transpose()).trace())
         ;
      TRACE_STRUCT.traceconjTYeTpYe = Re((TYe.conjugate()*Ye.transpose()).trace())
         ;
      TRACE_STRUCT.traceconjTYuTpYu = Re((TYu.conjugate()*Yu.transpose()).trace())
         ;
      TRACE_STRUCT.traceconjTYvTpYv = Re((TYv.conjugate()*Yv.transpose()).trace())
         ;
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYdAdjYdTYdAdjYd = Re((Yd*Yd.adjoint()*TYd*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd = Re((Yd*Yd.adjoint()*TYd*(TYd).adjoint()
         ).trace());
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYdAdjYuTYuAdjYd = Re((Yd*Yu.adjoint()*TYu*Yd.adjoint())
         .trace());
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd = Re((Yd*Yu.adjoint()*TYu*(TYd).adjoint()
         ).trace());
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd = Re((Yd*(TYd).adjoint()*TYd*Yd.adjoint()
         ).trace());
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd = Re((Yd*(TYu).adjoint()*TYu*Yd.adjoint()
         ).trace());
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
         .trace());
      TRACE_STRUCT.traceYeAdjYeTYeAdjYe = Re((Ye*Ye.adjoint()*TYe*Ye.adjoint())
         .trace());
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe = Re((Ye*Ye.adjoint()*TYe*(TYe).adjoint()
         ).trace());
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe = Re((Ye*(TYe).adjoint()*TYe*Ye.adjoint()
         ).trace());
      TRACE_STRUCT.traceYeconjTYvTpTYvAdjYe = Re((Ye*TYv.conjugate()*(TYv)
         .transpose()*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYdTYdAdjYu = Re((Yu*Yd.adjoint()*TYd*Yu.adjoint())
         .trace());
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu = Re((Yu*Yd.adjoint()*TYd*(TYu).adjoint()
         ).trace());
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
         .trace());
      TRACE_STRUCT.traceYuAdjYuTYuAdjYu = Re((Yu*Yu.adjoint()*TYu*Yu.adjoint())
         .trace());
      TRACE_STRUCT.traceYuAdjYuTYuAdjTYu = Re((Yu*Yu.adjoint()*TYu*(TYu).adjoint()
         ).trace());
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu = Re((Yu*(TYd).adjoint()*TYd*Yu.adjoint()
         ).trace());
      TRACE_STRUCT.traceYuAdjTYuTYuAdjYu = Re((Yu*(TYu).adjoint()*TYu*Yu.adjoint()
         ).trace());
      TRACE_STRUCT.traceYvAdjYvYvAdjYv = Re((Yv*Yv.adjoint()*Yv*Yv.adjoint())
         .trace());
      TRACE_STRUCT.traceYvAdjYvTYvAdjYv = Re((Yv*Yv.adjoint()*TYv*Yv.adjoint())
         .trace());
      TRACE_STRUCT.traceYvAdjYvTYvAdjTYv = Re((Yv*Yv.adjoint()*TYv*(TYv).adjoint()
         ).trace());
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe = Re((Yv*Yv.adjoint()*Ye.transpose()*
         Ye.conjugate()).trace());
      TRACE_STRUCT.traceYvAdjYvTpTYeconjTYe = Re((Yv*Yv.adjoint()*(TYe).transpose(
         )*TYe.conjugate()).trace());
      TRACE_STRUCT.traceYvAdjTYvTYvAdjYv = Re((Yv*(TYv).adjoint()*TYv*Yv.adjoint()
         ).trace());
      TRACE_STRUCT.traceAdjYeTYeconjYvTpYv = Re((Ye.adjoint()*TYe*Yv.conjugate()*
         Yv.transpose()).trace());
      TRACE_STRUCT.traceAdjYeTYeconjTYvTpYv = Re((Ye.adjoint()*TYe*TYv.conjugate()
         *Yv.transpose()).trace());
      TRACE_STRUCT.traceAdjYvTpYeconjYeTYv = Re((Yv.adjoint()*Ye.transpose()*
         Ye.conjugate()*TYv).trace());
      TRACE_STRUCT.traceAdjYvTpYeconjTYeTYv = Re((Yv.adjoint()*Ye.transpose()*
         TYe.conjugate()*TYv).trace());
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd = Re((md2*Yd*Yd.adjoint()*Yd*Yd.adjoint(
         )).trace());
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd = Re((md2*Yd*Yu.adjoint()*Yu*Yd.adjoint(
         )).trace());
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe = Re((me2*Ye*Ye.adjoint()*Ye*Ye.adjoint(
         )).trace());
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe = Re((ml2*Ye.adjoint()*Ye*Ye.adjoint()*
         Ye).trace());
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd = Re((mq2*Yd.adjoint()*Yd*Yd.adjoint()*
         Yd).trace());
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu = Re((mq2*Yd.adjoint()*Yd*Yu.adjoint()*
         Yu).trace());
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd = Re((mq2*Yu.adjoint()*Yu*Yd.adjoint()*
         Yd).trace());
      TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu = Re((mq2*Yu.adjoint()*Yu*Yu.adjoint()*
         Yu).trace());
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu = Re((mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint(
         )).trace());
      TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu = Re((mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint(
         )).trace());
      TRACE_STRUCT.traceYvAdjYvYvAdjYvconjml2 = Re((Yv*Yv.adjoint()*Yv*Yv.adjoint(
         )*ml2.conjugate()).trace());
      TRACE_STRUCT.traceYvAdjYvYvconjmvR2AdjYv = Re((Yv*Yv.adjoint()*Yv*
         mvR2.conjugate()*Yv.adjoint()).trace());
      TRACE_STRUCT.traceYvAdjYvconjml2YvAdjYv = Re((Yv*Yv.adjoint()*ml2.conjugate(
         )*Yv*Yv.adjoint()).trace());
      TRACE_STRUCT.traceYvAdjYvconjml2TpYeconjYe = Re((Yv*Yv.adjoint()*
         ml2.conjugate()*Ye.transpose()*Ye.conjugate()).trace());
      TRACE_STRUCT.traceYvAdjYvTpYeconjme2conjYe = Re((Yv*Yv.adjoint()*
         Ye.transpose()*me2.conjugate()*Ye.conjugate()).trace());
      TRACE_STRUCT.traceYvAdjYvTpYeconjYeconjml2 = Re((Yv*Yv.adjoint()*
         Ye.transpose()*Ye.conjugate()*ml2.conjugate()).trace());
      TRACE_STRUCT.traceYvconjmvR2AdjYvYvAdjYv = Re((Yv*mvR2.conjugate()*
         Yv.adjoint()*Yv*Yv.adjoint()).trace());
      TRACE_STRUCT.traceYvconjmvR2AdjYvTpYeconjYe = Re((Yv*mvR2.conjugate()*
         Yv.adjoint()*Ye.transpose()*Ye.conjugate()).trace());

   }

   if (get_loops() > 2) {

   }
}

std::ostream& operator<<(std::ostream& ostr, const UMSSM_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
