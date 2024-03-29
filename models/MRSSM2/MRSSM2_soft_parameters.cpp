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


#include "MRSSM2_soft_parameters.hpp"
#include "wrappers.hpp"
#include "functors.hpp"

#include <iostream>

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces
#define TRACE_STRUCT_TYPE Soft_traces
#define CALCULATE_TRACES(l) calc_soft_traces(l);

const int MRSSM2_soft_parameters::numberOfParameters;

MRSSM2_soft_parameters::MRSSM2_soft_parameters(const MRSSM2_input_parameters& input_)
   : MRSSM2_susy_parameters(input_)
{
   set_number_of_parameters(numberOfParameters);
}

MRSSM2_soft_parameters::MRSSM2_soft_parameters(
   const MRSSM2_susy_parameters& susy_model
   , double BMu_, double BMuD_, double BMuU_, const Eigen::Matrix<double,3,3>&
   mq2_, const Eigen::Matrix<double,3,3>& ml2_, double mHd2_, double mHu2_,
   const Eigen::Matrix<double,3,3>& md2_, const Eigen::Matrix<double,3,3>& mu2_
   , const Eigen::Matrix<double,3,3>& me2_, double mS2_, double mT2_, double
   moc2_, double mRd2_, double mRu2_, double MDBS_, double MDWBT_, double
   MDGoc_
)
   : MRSSM2_susy_parameters(susy_model)
   , BMu(BMu_), BMuD(BMuD_), BMuU(BMuU_), mq2(mq2_), ml2(ml2_), mHd2(mHd2_), mHu2(
   mHu2_), md2(md2_), mu2(mu2_), me2(me2_), mS2(mS2_), mT2(mT2_), moc2(moc2_),
   mRd2(mRd2_), mRu2(mRu2_), MDBS(MDBS_), MDWBT(MDWBT_), MDGoc(MDGoc_)
{
   set_number_of_parameters(numberOfParameters);
}

Eigen::ArrayXd MRSSM2_soft_parameters::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

MRSSM2_soft_parameters MRSSM2_soft_parameters::calc_beta(int loops) const
{
   double beta_BMu = 0.;
   double beta_BMuD = 0.;
   double beta_BMuU = 0.;
   Eigen::Matrix<double,3,3> beta_mq2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_ml2 = Eigen::Matrix<double,3,3>::Zero();
   double beta_mHd2 = 0.;
   double beta_mHu2 = 0.;
   Eigen::Matrix<double,3,3> beta_md2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_mu2 = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_me2 = Eigen::Matrix<double,3,3>::Zero();
   double beta_mS2 = 0.;
   double beta_mT2 = 0.;
   double beta_moc2 = 0.;
   double beta_mRd2 = 0.;
   double beta_mRu2 = 0.;
   double beta_MDBS = 0.;
   double beta_MDWBT = 0.;
   double beta_MDGoc = 0.;

   if (loops > 0) {
      const auto TRACE_STRUCT = CALCULATE_TRACES(loops);

      beta_BMu += calc_beta_BMu_1_loop(TRACE_STRUCT);
      beta_BMuD += calc_beta_BMuD_1_loop(TRACE_STRUCT);
      beta_BMuU += calc_beta_BMuU_1_loop(TRACE_STRUCT);
      beta_mq2 += calc_beta_mq2_1_loop(TRACE_STRUCT);
      beta_ml2 += calc_beta_ml2_1_loop(TRACE_STRUCT);
      beta_mHd2 += calc_beta_mHd2_1_loop(TRACE_STRUCT);
      beta_mHu2 += calc_beta_mHu2_1_loop(TRACE_STRUCT);
      beta_md2 += calc_beta_md2_1_loop(TRACE_STRUCT);
      beta_mu2 += calc_beta_mu2_1_loop(TRACE_STRUCT);
      beta_me2 += calc_beta_me2_1_loop(TRACE_STRUCT);
      beta_mS2 += calc_beta_mS2_1_loop(TRACE_STRUCT);
      beta_mT2 += calc_beta_mT2_1_loop(TRACE_STRUCT);
      beta_moc2 += calc_beta_moc2_1_loop(TRACE_STRUCT);
      beta_mRd2 += calc_beta_mRd2_1_loop(TRACE_STRUCT);
      beta_mRu2 += calc_beta_mRu2_1_loop(TRACE_STRUCT);
      beta_MDBS += calc_beta_MDBS_1_loop(TRACE_STRUCT);
      beta_MDWBT += calc_beta_MDWBT_1_loop(TRACE_STRUCT);
      beta_MDGoc += calc_beta_MDGoc_1_loop(TRACE_STRUCT);

      if (loops > 1) {
         beta_BMu += calc_beta_BMu_2_loop(TRACE_STRUCT);
         beta_BMuD += calc_beta_BMuD_2_loop(TRACE_STRUCT);
         beta_BMuU += calc_beta_BMuU_2_loop(TRACE_STRUCT);
         beta_mq2 += calc_beta_mq2_2_loop(TRACE_STRUCT);
         beta_ml2 += calc_beta_ml2_2_loop(TRACE_STRUCT);
         beta_mHd2 += calc_beta_mHd2_2_loop(TRACE_STRUCT);
         beta_mHu2 += calc_beta_mHu2_2_loop(TRACE_STRUCT);
         beta_md2 += calc_beta_md2_2_loop(TRACE_STRUCT);
         beta_mu2 += calc_beta_mu2_2_loop(TRACE_STRUCT);
         beta_me2 += calc_beta_me2_2_loop(TRACE_STRUCT);
         beta_mS2 += calc_beta_mS2_2_loop(TRACE_STRUCT);
         beta_mT2 += calc_beta_mT2_2_loop(TRACE_STRUCT);
         beta_moc2 += calc_beta_moc2_2_loop(TRACE_STRUCT);
         beta_mRd2 += calc_beta_mRd2_2_loop(TRACE_STRUCT);
         beta_mRu2 += calc_beta_mRu2_2_loop(TRACE_STRUCT);
         beta_MDBS += calc_beta_MDBS_2_loop(TRACE_STRUCT);
         beta_MDWBT += calc_beta_MDWBT_2_loop(TRACE_STRUCT);
         beta_MDGoc += calc_beta_MDGoc_2_loop(TRACE_STRUCT);

         if (loops > 2) {

            if (loops > 3) {

               if (loops > 4) {

               }
            }
         }
      }
   }


   const MRSSM2_susy_parameters susy_betas(MRSSM2_susy_parameters::calc_beta(loops));

   return MRSSM2_soft_parameters(susy_betas, beta_BMu, beta_BMuD, beta_BMuU, beta_mq2, beta_ml2, beta_mHd2, beta_mHu2, beta_md2, beta_mu2, beta_me2, beta_mS2, beta_mT2, beta_moc2, beta_mRd2, beta_mRu2, beta_MDBS, beta_MDWBT, beta_MDGoc);
}

MRSSM2_soft_parameters MRSSM2_soft_parameters::calc_beta() const
{
   return calc_beta(get_loops());
}

void MRSSM2_soft_parameters::clear()
{
   MRSSM2_susy_parameters::clear();

   BMu = 0.;
   BMuD = 0.;
   BMuU = 0.;
   mq2 = Eigen::Matrix<double,3,3>::Zero();
   ml2 = Eigen::Matrix<double,3,3>::Zero();
   mHd2 = 0.;
   mHu2 = 0.;
   md2 = Eigen::Matrix<double,3,3>::Zero();
   mu2 = Eigen::Matrix<double,3,3>::Zero();
   me2 = Eigen::Matrix<double,3,3>::Zero();
   mS2 = 0.;
   mT2 = 0.;
   moc2 = 0.;
   mRd2 = 0.;
   mRu2 = 0.;
   MDBS = 0.;
   MDWBT = 0.;
   MDGoc = 0.;

}

Eigen::ArrayXd MRSSM2_soft_parameters::get() const
{
   Eigen::ArrayXd pars(MRSSM2_susy_parameters::get());
   pars.conservativeResize(numberOfParameters);

   pars(41) = BMu;
   pars(42) = BMuD;
   pars(43) = BMuU;
   pars(44) = mq2(0,0);
   pars(45) = mq2(0,1);
   pars(46) = mq2(0,2);
   pars(47) = mq2(1,0);
   pars(48) = mq2(1,1);
   pars(49) = mq2(1,2);
   pars(50) = mq2(2,0);
   pars(51) = mq2(2,1);
   pars(52) = mq2(2,2);
   pars(53) = ml2(0,0);
   pars(54) = ml2(0,1);
   pars(55) = ml2(0,2);
   pars(56) = ml2(1,0);
   pars(57) = ml2(1,1);
   pars(58) = ml2(1,2);
   pars(59) = ml2(2,0);
   pars(60) = ml2(2,1);
   pars(61) = ml2(2,2);
   pars(62) = mHd2;
   pars(63) = mHu2;
   pars(64) = md2(0,0);
   pars(65) = md2(0,1);
   pars(66) = md2(0,2);
   pars(67) = md2(1,0);
   pars(68) = md2(1,1);
   pars(69) = md2(1,2);
   pars(70) = md2(2,0);
   pars(71) = md2(2,1);
   pars(72) = md2(2,2);
   pars(73) = mu2(0,0);
   pars(74) = mu2(0,1);
   pars(75) = mu2(0,2);
   pars(76) = mu2(1,0);
   pars(77) = mu2(1,1);
   pars(78) = mu2(1,2);
   pars(79) = mu2(2,0);
   pars(80) = mu2(2,1);
   pars(81) = mu2(2,2);
   pars(82) = me2(0,0);
   pars(83) = me2(0,1);
   pars(84) = me2(0,2);
   pars(85) = me2(1,0);
   pars(86) = me2(1,1);
   pars(87) = me2(1,2);
   pars(88) = me2(2,0);
   pars(89) = me2(2,1);
   pars(90) = me2(2,2);
   pars(91) = mS2;
   pars(92) = mT2;
   pars(93) = moc2;
   pars(94) = mRd2;
   pars(95) = mRu2;
   pars(96) = MDBS;
   pars(97) = MDWBT;
   pars(98) = MDGoc;


   return pars;
}

void MRSSM2_soft_parameters::print(std::ostream& ostr) const
{
   MRSSM2_susy_parameters::print(ostr);
   ostr << "soft parameters at Q = " << get_scale() << ":\n";
   ostr << "BMu = " << BMu << '\n';
   ostr << "BMuD = " << BMuD << '\n';
   ostr << "BMuU = " << BMuU << '\n';
   ostr << "mq2 = " << mq2 << '\n';
   ostr << "ml2 = " << ml2 << '\n';
   ostr << "mHd2 = " << mHd2 << '\n';
   ostr << "mHu2 = " << mHu2 << '\n';
   ostr << "md2 = " << md2 << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "me2 = " << me2 << '\n';
   ostr << "mS2 = " << mS2 << '\n';
   ostr << "mT2 = " << mT2 << '\n';
   ostr << "moc2 = " << moc2 << '\n';
   ostr << "mRd2 = " << mRd2 << '\n';
   ostr << "mRu2 = " << mRu2 << '\n';
   ostr << "MDBS = " << MDBS << '\n';
   ostr << "MDWBT = " << MDWBT << '\n';
   ostr << "MDGoc = " << MDGoc << '\n';

}

void MRSSM2_soft_parameters::set(const Eigen::ArrayXd& pars)
{
   MRSSM2_susy_parameters::set(pars);

   BMu = pars(41);
   BMuD = pars(42);
   BMuU = pars(43);
   mq2(0,0) = pars(44);
   mq2(0,1) = pars(45);
   mq2(0,2) = pars(46);
   mq2(1,0) = pars(47);
   mq2(1,1) = pars(48);
   mq2(1,2) = pars(49);
   mq2(2,0) = pars(50);
   mq2(2,1) = pars(51);
   mq2(2,2) = pars(52);
   ml2(0,0) = pars(53);
   ml2(0,1) = pars(54);
   ml2(0,2) = pars(55);
   ml2(1,0) = pars(56);
   ml2(1,1) = pars(57);
   ml2(1,2) = pars(58);
   ml2(2,0) = pars(59);
   ml2(2,1) = pars(60);
   ml2(2,2) = pars(61);
   mHd2 = pars(62);
   mHu2 = pars(63);
   md2(0,0) = pars(64);
   md2(0,1) = pars(65);
   md2(0,2) = pars(66);
   md2(1,0) = pars(67);
   md2(1,1) = pars(68);
   md2(1,2) = pars(69);
   md2(2,0) = pars(70);
   md2(2,1) = pars(71);
   md2(2,2) = pars(72);
   mu2(0,0) = pars(73);
   mu2(0,1) = pars(74);
   mu2(0,2) = pars(75);
   mu2(1,0) = pars(76);
   mu2(1,1) = pars(77);
   mu2(1,2) = pars(78);
   mu2(2,0) = pars(79);
   mu2(2,1) = pars(80);
   mu2(2,2) = pars(81);
   me2(0,0) = pars(82);
   me2(0,1) = pars(83);
   me2(0,2) = pars(84);
   me2(1,0) = pars(85);
   me2(1,1) = pars(86);
   me2(1,2) = pars(87);
   me2(2,0) = pars(88);
   me2(2,1) = pars(89);
   me2(2,2) = pars(90);
   mS2 = pars(91);
   mT2 = pars(92);
   moc2 = pars(93);
   mRd2 = pars(94);
   mRu2 = pars(95);
   MDBS = pars(96);
   MDWBT = pars(97);
   MDGoc = pars(98);

}

MRSSM2_soft_parameters::Soft_traces MRSSM2_soft_parameters::calc_soft_traces(int loops) const
{
   Soft_traces soft_traces;

   if (loops > 0) {
      
      TRACE_STRUCT.Tr11 = Re(0.7745966692414834*g1*(-mHd2 + mHu2 + mRd2 - mRu2 + (md2
         ).trace() + (me2).trace() - (ml2).trace() + (mq2).trace() - 2*(mu2).trace())
         );
      TRACE_STRUCT.Tr2U111 = Re(0.1*Sqr(g1)*(3*mHd2 + 3*mHu2 + 3*mRd2 + 3*mRu2 + 2*(
         md2).trace() + 6*(me2).trace() + 3*(ml2).trace() + (mq2).trace() + 8*(mu2).
         trace()));
      TRACE_STRUCT.Tr31 = Re(0.012909944487358056*g1*(30*(mHd2 - mRd2)*AbsSqr(LamSD)
         - 30*(mHu2 - mRu2)*AbsSqr(LamSU) + 45*mHd2*AbsSqr(LamTD) - 45*mRd2*AbsSqr(
         LamTD) - 45*mHu2*AbsSqr(LamTU) + 45*mRu2*AbsSqr(LamTU) - 9*mHd2*Sqr(g1) + 9*
         mHu2*Sqr(g1) + 9*mRd2*Sqr(g1) - 9*mRu2*Sqr(g1) - 45*mHd2*Sqr(g2) + 45*mHu2*
         Sqr(g2) + 45*mRd2*Sqr(g2) - 45*mRu2*Sqr(g2) + 4*Sqr(g1)*(md2).trace() + 80*
         Sqr(g3)*(md2).trace() + 36*Sqr(g1)*(me2).trace() - 9*Sqr(g1)*(ml2).trace() -
         45*Sqr(g2)*(ml2).trace() + Sqr(g1)*(mq2).trace() + 45*Sqr(g2)*(mq2).trace()
         + 80*Sqr(g3)*(mq2).trace() - 32*Sqr(g1)*(mu2).trace() - 160*Sqr(g3)*(mu2).
         trace() + 90*mHd2*(Yd*Yd.adjoint()).trace() + 30*mHd2*(Ye*Ye.adjoint()).
         trace() - 90*mHu2*(Yu*Yu.adjoint()).trace() - 60*(Yd*Yd.adjoint()*md2.
         conjugate()).trace() - 30*(Yd*mq2.conjugate()*Yd.adjoint()).trace() - 60*(Ye
         *Ye.adjoint()*me2.conjugate()).trace() + 30*(Ye*ml2.conjugate()*Ye.adjoint()
         ).trace() + 120*(Yu*Yu.adjoint()*mu2.conjugate()).trace() - 30*(Yu*mq2.
         conjugate()*Yu.adjoint()).trace()));
      TRACE_STRUCT.Tr22 = Re(0.5*(mHd2 + mHu2 + mRd2 + mRu2 + 4*mT2 + (ml2).trace() +
         3*(mq2).trace()));
      TRACE_STRUCT.Tr23 = Re(0.5*(6*moc2 + (md2).trace() + 2*(mq2).trace() + (mu2).
         trace()));

      TRACE_STRUCT.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      TRACE_STRUCT.tracemd2YdAdjYd = Re((md2*Yd*Yd.adjoint()).trace());
      TRACE_STRUCT.traceme2YeAdjYe = Re((me2*Ye*Ye.adjoint()).trace());
      TRACE_STRUCT.traceml2AdjYeYe = Re((ml2*Ye.adjoint()*Ye).trace());
      TRACE_STRUCT.tracemq2AdjYdYd = Re((mq2*Yd.adjoint()*Yd).trace());
      TRACE_STRUCT.tracemq2AdjYuYu = Re((mq2*Yu.adjoint()*Yu).trace());
      TRACE_STRUCT.tracemu2YuAdjYu = Re((mu2*Yu*Yu.adjoint()).trace());

   }

   if (loops > 1) {
      TRACE_STRUCT.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint()).trace()
         );
      TRACE_STRUCT.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint()).trace()
         );
      TRACE_STRUCT.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint()).trace()
         );
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

std::ostream& operator<<(std::ostream& ostr, const MRSSM2_soft_parameters& soft_pars)
{
   soft_pars.print(ostr);
   return ostr;
}

} // namespace flexiblesusy
