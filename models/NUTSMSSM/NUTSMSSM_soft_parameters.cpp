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


#include "NUTSMSSM_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES(l) calc_soft_traces(l);

const int NUTSMSSM_soft_parameters::numberOfParameters;

NUTSMSSM_soft_parameters::NUTSMSSM_soft_parameters(const NUTSMSSM_input_parameters& input_)
   : NUTSMSSM_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

NUTSMSSM_soft_parameters::NUTSMSSM_soft_parameters(
   const NUTSMSSM_susy_parameters& susy_model
   , const Eigen::Matrix<double,3,3>& TYd_, const Eigen::Matrix<double,3,3>& TYe_,
   double TLambdax_, double TKappa_, const Eigen::Matrix<double,3,3>& TYu_,
   double BMu_, double BMS_, double LL1_, const Eigen::Matrix<double,3,3>& mq2_
   , const Eigen::Matrix<double,3,3>& ml2_, double mHd2_, double mHu2_, const
   Eigen::Matrix<double,3,3>& md2_, const Eigen::Matrix<double,3,3>& mu2_,
   const Eigen::Matrix<double,3,3>& me2_, double ms2_, double MassB_, double
   MassWB_, double MassG_
)
   : NUTSMSSM_susy_parameters(susy_model)
   , TYd(TYd_), TYe(TYe_), TLambdax(TLambdax_), TKappa(TKappa_), TYu(TYu_), BMu(
   BMu_), BMS(BMS_), LL1(LL1_), mq2(mq2_), ml2(ml2_), mHd2(mHd2_), mHu2(mHu2_),
   md2(md2_), mu2(mu2_), me2(me2_), ms2(ms2_), MassB(MassB_), MassWB(MassWB_),
   MassG(MassG_)
{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd NUTSMSSM_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

NUTSMSSM_soft_parameters NUTSMSSM_soft_parameters::calc_beta(int loops) const
{
   Eigen::Matrix<double,3,3> beta_TYd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_TYe = Eigen::Matrix<double,3,3>::Zero();
   double beta_TLambdax = 0.;
   double beta_TKappa = 0.;
   Eigen::Matrix<double,3,3> beta_TYu = Eigen::Matrix<double,3,3>::Zero();
   double beta_BMu = 0.;
   double beta_BMS = 0.;
   double beta_LL1 = 0.;
   Eigen::Matrix<double,3,3> beta_mq2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_ml2 = Eigen::Matrix<double,3,3>::Zero();
   double beta_mHd2 = 0.;
   double beta_mHu2 = 0.;
   Eigen::Matrix<double,3,3> beta_md2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_mu2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_me2 = Eigen::Matrix<double,3,3>::Zero();
   double beta_ms2 = 0.;
   double beta_MassB = 0.;
   double beta_MassWB = 0.;
   double beta_MassG = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_TYd += calc_beta_TYd_1_loop(TRACE_STRUCT);
      beta_TYe += calc_beta_TYe_1_loop(TRACE_STRUCT);
      beta_TLambdax += calc_beta_TLambdax_1_loop(TRACE_STRUCT);
      beta_TKappa += calc_beta_TKappa_1_loop(TRACE_STRUCT);
      beta_TYu += calc_beta_TYu_1_loop(TRACE_STRUCT);
      beta_BMu += calc_beta_BMu_1_loop(TRACE_STRUCT);
      beta_BMS += calc_beta_BMS_1_loop(TRACE_STRUCT);
      beta_LL1 += calc_beta_LL1_1_loop(TRACE_STRUCT);
      beta_mq2 += calc_beta_mq2_1_loop(TRACE_STRUCT);
      beta_ml2 += calc_beta_ml2_1_loop(TRACE_STRUCT);
      beta_mHd2 += calc_beta_mHd2_1_loop(TRACE_STRUCT);
      beta_mHu2 += calc_beta_mHu2_1_loop(TRACE_STRUCT);
      beta_md2 += calc_beta_md2_1_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_1_loop(TRACE_STRUCT);
      beta_me2 += calc_beta_me2_1_loop(TRACE_STRUCT);
      beta_ms2 += calc_beta_ms2_1_loop(TRACE_STRUCT);
      beta_MassB += calc_beta_MassB_1_loop(TRACE_STRUCT);
      beta_MassWB += calc_beta_MassWB_1_loop(TRACE_STRUCT);
      beta_MassG += calc_beta_MassG_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_TYd += calc_beta_TYd_2_loop(TRACE_STRUCT);
         beta_TYe += calc_beta_TYe_2_loop(TRACE_STRUCT);
         beta_TLambdax += calc_beta_TLambdax_2_loop(TRACE_STRUCT);
         beta_TKappa += calc_beta_TKappa_2_loop(TRACE_STRUCT);
         beta_TYu += calc_beta_TYu_2_loop(TRACE_STRUCT);
         beta_BMu += calc_beta_BMu_2_loop(TRACE_STRUCT);
         beta_BMS += calc_beta_BMS_2_loop(TRACE_STRUCT);
         beta_LL1 += calc_beta_LL1_2_loop(TRACE_STRUCT);
         beta_mq2 += calc_beta_mq2_2_loop(TRACE_STRUCT);
         beta_ml2 += calc_beta_ml2_2_loop(TRACE_STRUCT);
         beta_mHd2 += calc_beta_mHd2_2_loop(TRACE_STRUCT);
         beta_mHu2 += calc_beta_mHu2_2_loop(TRACE_STRUCT);
         beta_md2 += calc_beta_md2_2_loop(TRACE_STRUCT);
         beta_mu2 += calc_beta_mu2_2_loop(TRACE_STRUCT);
         beta_me2 += calc_beta_me2_2_loop(TRACE_STRUCT);
         beta_ms2 += calc_beta_ms2_2_loop(TRACE_STRUCT);
         beta_MassB += calc_beta_MassB_2_loop(TRACE_STRUCT);
         beta_MassWB += calc_beta_MassWB_2_loop(TRACE_STRUCT);
         beta_MassG += calc_beta_MassG_2_loop(TRACE_STRUCT);

         if (loops > 2) {

            if (loops > 3) {

               if (loops > 4) {

               }
            }
         }
      }
   }


   const NUTSMSSM_susy_parameters susy_betas(NUTSMSSM_susy_parameters::calc_beta(loops));

   return NUTSMSSM_soft_parameters(susy_betas, beta_TYd, beta_TYe, beta_TLambdax, beta_TKappa, beta_TYu, beta_BMu, beta_BMS, beta_LL1, beta_mq2, beta_ml2, beta_mHd2, beta_mHu2, beta_md2, beta_mu2, beta_me2, beta_ms2, beta_MassB, beta_MassWB, beta_MassG);
}

NUTSMSSM_soft_parameters NUTSMSSM_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void NUTSMSSM_soft_parameters::clear()
{
   NUTSMSSM_susy_parameters::clear();

   TYd = Eigen::Matrix<double,3,3>::Zero();
   TYe = Eigen::Matrix<double,3,3>::Zero();
   TLambdax = 0.;
   TKappa = 0.;
   TYu = Eigen::Matrix<double,3,3>::Zero();
   BMu = 0.;
   BMS = 0.;
   LL1 = 0.;
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   ms2 = 0.;
   MassB = 0.;
   MassWB = 0.;
   MassG = 0.;

}

Eigen::ArrayXd NUTSMSSM_soft_parameters::get() const
{
   Eigen::ArrayXd pars(NUTSMSSM_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(38) = TYd(0,0);
   pars(39) = TYd(0,1);
   pars(40) = TYd(0,2);
   pars(41) = TYd(1,0);
   pars(42) = TYd(1,1);
   pars(43) = TYd(1,2);
   pars(44) = TYd(2,0);
   pars(45) = TYd(2,1);
   pars(46) = TYd(2,2);
   pars(47) = TYe(0,0);
   pars(48) = TYe(0,1);
   pars(49) = TYe(0,2);
   pars(50) = TYe(1,0);
   pars(51) = TYe(1,1);
   pars(52) = TYe(1,2);
   pars(53) = TYe(2,0);
   pars(54) = TYe(2,1);
   pars(55) = TYe(2,2);
   pars(56) = TLambdax;
   pars(57) = TKappa;
   pars(58) = TYu(0,0);
   pars(59) = TYu(0,1);
   pars(60) = TYu(0,2);
   pars(61) = TYu(1,0);
   pars(62) = TYu(1,1);
   pars(63) = TYu(1,2);
   pars(64) = TYu(2,0);
   pars(65) = TYu(2,1);
   pars(66) = TYu(2,2);
   pars(67) = BMu;
   pars(68) = BMS;
   pars(69) = LL1;
   pars(70) = mq2(0,0);
   pars(71) = mq2(0,1);
   pars(72) = mq2(0,2);
   pars(73) = mq2(1,0);
   pars(74) = mq2(1,1);
   pars(75) = mq2(1,2);
   pars(76) = mq2(2,0);
   pars(77) = mq2(2,1);
   pars(78) = mq2(2,2);
   pars(79) = ml2(0,0);
   pars(80) = ml2(0,1);
   pars(81) = ml2(0,2);
   pars(82) = ml2(1,0);
   pars(83) = ml2(1,1);
   pars(84) = ml2(1,2);
   pars(85) = ml2(2,0);
   pars(86) = ml2(2,1);
   pars(87) = ml2(2,2);
   pars(88) = mHd2;
   pars(89) = mHu2;
   pars(90) = md2(0,0);
   pars(91) = md2(0,1);
   pars(92) = md2(0,2);
   pars(93) = md2(1,0);
   pars(94) = md2(1,1);
   pars(95) = md2(1,2);
   pars(96) = md2(2,0);
   pars(97) = md2(2,1);
   pars(98) = md2(2,2);
   pars(99) = mu2(0,0);
   pars(100) = mu2(0,1);
   pars(101) = mu2(0,2);
   pars(102) = mu2(1,0);
   pars(103) = mu2(1,1);
   pars(104) = mu2(1,2);
   pars(105) = mu2(2,0);
   pars(106) = mu2(2,1);
   pars(107) = mu2(2,2);
   pars(108) = me2(0,0);
   pars(109) = me2(0,1);
   pars(110) = me2(0,2);
   pars(111) = me2(1,0);
   pars(112) = me2(1,1);
   pars(113) = me2(1,2);
   pars(114) = me2(2,0);
   pars(115) = me2(2,1);
   pars(116) = me2(2,2);
   pars(117) = ms2;
   pars(118) = MassB;
   pars(119) = MassWB;
   pars(120) = MassG;


   return pars;
}

void NUTSMSSM_soft_parameters::print(std::ostream& ostr) const
{
   NUTSMSSM_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "TYd = " << TYd << '\n';
   ostr << "TYe = " << TYe << '\n';
   ostr << "TLambdax = " << TLambdax << '\n';
   ostr << "TKappa = " << TKappa << '\n';
   ostr << "TYu = " << TYu << '\n';
   ostr << "BMu = " << BMu << '\n';
   ostr << "BMS = " << BMS << '\n';
   ostr << "LL1 = " << LL1 << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "ms2 = " << ms2 << '\n';
   ostr << "MassB = " << MassB << '\n';
   ostr << "MassWB = " << MassWB << '\n';
   ostr << "MassG = " << MassG << '\n';

}

void NUTSMSSM_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   NUTSMSSM_susy_parameters::set(pars);

   TYd(0,0) = pars(38);
   TYd(0,1) = pars(39);
   TYd(0,2) = pars(40);
   TYd(1,0) = pars(41);
   TYd(1,1) = pars(42);
   TYd(1,2) = pars(43);
   TYd(2,0) = pars(44);
   TYd(2,1) = pars(45);
   TYd(2,2) = pars(46);
   TYe(0,0) = pars(47);
   TYe(0,1) = pars(48);
   TYe(0,2) = pars(49);
   TYe(1,0) = pars(50);
   TYe(1,1) = pars(51);
   TYe(1,2) = pars(52);
   TYe(2,0) = pars(53);
   TYe(2,1) = pars(54);
   TYe(2,2) = pars(55);
   TLambdax = pars(56);
   TKappa = pars(57);
   TYu(0,0) = pars(58);
   TYu(0,1) = pars(59);
   TYu(0,2) = pars(60);
   TYu(1,0) = pars(61);
   TYu(1,1) = pars(62);
   TYu(1,2) = pars(63);
   TYu(2,0) = pars(64);
   TYu(2,1) = pars(65);
   TYu(2,2) = pars(66);
   BMu = pars(67);
   BMS = pars(68);
   LL1 = pars(69);
   mq2(0,0) = pars(70);
   mq2(0,1) = pars(71);
   mq2(0,2) = pars(72);
   mq2(1,0) = pars(73);
   mq2(1,1) = pars(74);
   mq2(1,2) = pars(75);
   mq2(2,0) = pars(76);
   mq2(2,1) = pars(77);
   mq2(2,2) = pars(78);
   ml2(0,0) = pars(79);
   ml2(0,1) = pars(80);
   ml2(0,2) = pars(81);
   ml2(1,0) = pars(82);
   ml2(1,1) = pars(83);
   ml2(1,2) = pars(84);
   ml2(2,0) = pars(85);
   ml2(2,1) = pars(86);
   ml2(2,2) = pars(87);
   mHd2 = pars(88);
   mHu2 = pars(89);
   md2(0,0) = pars(90);
   md2(0,1) = pars(91);
   md2(0,2) = pars(92);
   md2(1,0) = pars(93);
   md2(1,1) = pars(94);
   md2(1,2) = pars(95);
   md2(2,0) = pars(96);
   md2(2,1) = pars(97);
   md2(2,2) = pars(98);
   mu2(0,0) = pars(99);
   mu2(0,1) = pars(100);
   mu2(0,2) = pars(101);
   mu2(1,0) = pars(102);
   mu2(1,1) = pars(103);
   mu2(1,2) = pars(104);
   mu2(2,0) = pars(105);
   mu2(2,1) = pars(106);
   mu2(2,2) = pars(107);
   me2(0,0) = pars(108);
   me2(0,1) = pars(109);
   me2(0,2) = pars(110);
   me2(1,0) = pars(111);
   me2(1,1) = pars(112);
   me2(1,2) = pars(113);
   me2(2,0) = pars(114);
   me2(2,1) = pars(115);
   me2(2,2) = pars(116);
   ms2 = pars(117);
   MassB = pars(118);
   MassWB = pars(119);
   MassG = pars(120);

}

NUTSMSSM_soft_parameters::Soft_traces NUTSMSSM_soft_parameters::calc_soft_traces(int loops) const
{
   Soft_traces soft_traces;

   if (loops > 0) {
      
      TRACE_STRUCT.Tr11 = Re(0.7745966692414834*g1*(-mHd2 + mHu2 + (md2).trace() + (
         me2).trace() - (ml2).trace() + (mq2).trace() - 2*(mu2).trace()));
      TRACE_STRUCT.Tr2U111 = Re(0.1*Sqr(g1)*(3*mHd2 + 3*mHu2 + 2*(md2).trace() + 6*(
         me2).trace() + 3*(ml2).trace() + (mq2).trace() + 8*(mu2).trace()));
      TRACE_STRUCT.Tr31 = Re(0.012909944487358056*g1*(30*(mHd2 - mHu2)*AbsSqr(Lambdax
         ) - 9*mHd2*Sqr(g1) + 9*mHu2*Sqr(g1) - 45*mHd2*Sqr(g2) + 45*mHu2*Sqr(g2) + 4*
         (Sqr(g1) + 20*Sqr(g3))*(md2).trace() + 36*Sqr(g1)*(me2).trace() - 9*Sqr(g1)*
         (ml2).trace() - 45*Sqr(g2)*(ml2).trace() + Sqr(g1)*(mq2).trace() + 45*Sqr(g2
         )*(mq2).trace() + 80*Sqr(g3)*(mq2).trace() - 32*Sqr(g1)*(mu2).trace() - 160*
         Sqr(g3)*(mu2).trace() + 90*mHd2*(Yd*Yd.adjoint()).trace() + 30*mHd2*(Ye*Ye.
         adjoint()).trace() - 90*mHu2*(Yu*Yu.adjoint()).trace() - 60*(Yd*Yd.adjoint()
         *md2.conjugate()).trace() - 30*(Yd*mq2.conjugate()*Yd.adjoint()).trace() -
         60*(Ye*Ye.adjoint()*me2.conjugate()).trace() + 30*(Ye*ml2.conjugate()*Ye.
         adjoint()).trace() + 120*(Yu*Yu.adjoint()*mu2.conjugate()).trace() - 30*(Yu*
         mq2.conjugate()*Yu.adjoint()).trace()));
      TRACE_STRUCT.Tr22 = Re(0.5*(mHd2 + mHu2 + (ml2).trace() + 3*(mq2).trace()));
      TRACE_STRUCT.Tr23 = Re(0.5*((md2).trace() + 2*(mq2).trace() + (mu2).trace()));

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.traceAdjYdTYd = Re((Yd.adjoint()*TYd).trace());
      TRACE_STRUCT.traceAdjYeTYe = Re((Ye.adjoint()*TYe).trace());
      TRACE_STRUCT.traceAdjYuTYu = Re((Yu.adjoint()*TYu).trace());
      TRACE_STRUCT.traceconjTYdTpTYd = Re((TYd.conjugate()*(TYd).transpose()).trace()
         );
      TRACE_STRUCT.traceconjTYeTpTYe = Re((TYe.conjugate()*(TYe).transpose()).trace()
         );
      TRACE_STRUCT.traceconjTYuTpTYu = Re((TYu.conjugate()*(TYu).transpose()).trace()
         );
      TRACE_STRUCT.tracemd2YdAdjYd = Re((md2*Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceme2YeAdjYe = Re((me2*Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceml2AdjYeYe = Re((ml2*Ye.adjoint()*Ye).trace());
      TRACE_STRUCT.tracemq2AdjYdYd = Re((mq2*Yd.adjoint()*Yd).trace());
      TRACE_STRUCT.tracemq2AdjYuYu = Re((mq2*Yu.adjoint()*Yu).trace());
      TRACE_STRUCT.tracemu2YuAdjYu = Re((mu2*Yu*Yu.adjoint()).trace());

   }

   if (loops > 1) {
      TRACE_STRUCT.traceconjTYdTpYd = Re((TYd.conjugate()*Yd.transpose()).trace());
      TRACE_STRUCT.traceconjTYeTpYe = Re((TYe.conjugate()*Ye.transpose()).trace());
      TRACE_STRUCT.traceconjTYuTpYu = Re((TYu.conjugate()*Yu.transpose()).trace());
      TRACE_STRUCT.traceYdAdjYdconjmd2 = Re((Yd*Yd.adjoint()*md2.conjugate()).trace()
         );
      TRACE_STRUCT.traceYdconjmq2AdjYd = Re((Yd*mq2.conjugate()*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYeconjme2 = Re((Ye*Ye.adjoint()*me2.conjugate()).trace()
         );
      TRACE_STRUCT.traceYeconjml2AdjYe = Re((Ye*ml2.conjugate()*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYuAdjYuconjmu2 = Re((Yu*Yu.adjoint()*mu2.conjugate()).trace()
         );
      TRACE_STRUCT.traceYuconjmq2AdjYu = Re((Yu*mq2.conjugate()*Yu.adjoint()).trace()
         );
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
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd = Re((md2*Yd*Yd.adjoint()*Yd*Yd.adjoint()).
         trace());
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd = Re((md2*Yd*Yu.adjoint()*Yu*Yd.adjoint()).
         trace());
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe = Re((me2*Ye*Ye.adjoint()*Ye*Ye.adjoint()).
         trace());
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

   }

   if (loops > 2) {

   }

   return soft_traces;
}

std::ostream& operator<<(std::ostream& ostr, const NUTSMSSM_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
