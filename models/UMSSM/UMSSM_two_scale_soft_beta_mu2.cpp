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

// File generated at Sun 18 Oct 2015 12:18:56

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mu2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qu = INPUT(Qu);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (oneOver16PiSqr*(4*mHu2*(Yu*Yu.adjoint()) + 4*(TYu*(TYu)
      .adjoint()) + 2*(mu2*Yu*Yu.adjoint()) + 4*(Yu*mq2*Yu.adjoint()) + 2*(Yu*
      Yu.adjoint()*mu2) - 1.0327955589886444*g1*Tr11*UNITMATRIX(3) + 2*gp*Qu*
      Tr14*UNITMATRIX(3) - 2.1333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(
      3) - 10.666666666666666*AbsSqr(MassG)*Sqr(g3)*UNITMATRIX(3) - 8*AbsSqr(
      MassU)*Sqr(gp)*Sqr(Qu)*UNITMATRIX(3))).real();


   return beta_mu2;
}

/**
 * Calculates the two-loop beta function of mu2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (0.008888888888888889*twoLoop*(2*Conj(MassB)*Sqr(g1)*(45*(
      -2*MassB*(Yu*Yu.adjoint()) + TYu*Yu.adjoint()) + 4*(642*MassB*Sqr(g1) + 5
      *(16*(2*MassB + MassG)*Sqr(g3) - 3*(2*MassB + MassU)*(9*Qd + 9*Qe - 3*QHd
      + 3*QHu - 9*Ql + 9*Qq - 22*Qu)*Qu*Sqr(gp)))*UNITMATRIX(3)) + 5*(32*Conj(
      MassG)*Sqr(g3)*(4*(MassB + 2*MassG)*Sqr(g1) + 15*(-2*MassG*Sqr(g3) + (2*
      MassG + MassU)*Sqr(gp)*Sqr(Qu)))*UNITMATRIX(3) - 3*(3*(2*(15*
      traceconjTYuTpTYu + 15*tracemq2AdjYuYu + 15*tracemu2YuAdjYu + 30*mHu2*
      traceYuAdjYu + 5*(mHd2 + 2*mHu2 + ms2)*AbsSqr(Lambdax) + 5*AbsSqr(
      TLambdax) + mHu2*Sqr(g1) - 15*mHu2*Sqr(g2) - 30*AbsSqr(MassWB)*Sqr(g2) -
      10*mHu2*Sqr(gp)*Sqr(QHu) - 10*mHu2*Sqr(gp)*Sqr(Qq) + 10*mHu2*Sqr(gp)*Sqr(
      Qu))*(Yu*Yu.adjoint()) + (30*traceAdjYuTYu - 2*MassB*Sqr(g1) + 30*MassWB*
      Sqr(g2) + 20*MassU*Sqr(gp)*Sqr(QHu) + 20*MassU*Sqr(gp)*Sqr(Qq) - 20*MassU
      *Sqr(gp)*Sqr(Qu) + 10*Conj(Lambdax)*TLambdax)*(Yu*(TYu).adjoint()) + 30*
      traceconjTYuTpYu*(TYu*Yu.adjoint()) + 10*Conj(TLambdax)*Lambdax*(TYu*
      Yu.adjoint()) + 30*Conj(MassWB)*Sqr(g2)*(TYu*Yu.adjoint()) + 30*
      traceYuAdjYu*(TYu*(TYu).adjoint()) + 10*AbsSqr(Lambdax)*(TYu*(TYu)
      .adjoint()) + 2*Sqr(g1)*(TYu*(TYu).adjoint()) - 30*Sqr(g2)*(TYu*(TYu)
      .adjoint()) - 20*Sqr(gp)*Sqr(QHu)*(TYu*(TYu).adjoint()) - 20*Sqr(gp)*Sqr(
      Qq)*(TYu*(TYu).adjoint()) + 20*Sqr(gp)*Sqr(Qu)*(TYu*(TYu).adjoint()) + 15
      *traceYuAdjYu*(mu2*Yu*Yu.adjoint()) + 5*AbsSqr(Lambdax)*(mu2*Yu*
      Yu.adjoint()) + Sqr(g1)*(mu2*Yu*Yu.adjoint()) - 15*Sqr(g2)*(mu2*Yu*
      Yu.adjoint()) - 10*Sqr(gp)*Sqr(QHu)*(mu2*Yu*Yu.adjoint()) - 10*Sqr(gp)*
      Sqr(Qq)*(mu2*Yu*Yu.adjoint()) + 10*Sqr(gp)*Sqr(Qu)*(mu2*Yu*Yu.adjoint())
      + 30*traceYuAdjYu*(Yu*mq2*Yu.adjoint()) + 10*AbsSqr(Lambdax)*(Yu*mq2*
      Yu.adjoint()) + 2*Sqr(g1)*(Yu*mq2*Yu.adjoint()) - 30*Sqr(g2)*(Yu*mq2*
      Yu.adjoint()) - 20*Sqr(gp)*Sqr(QHu)*(Yu*mq2*Yu.adjoint()) - 20*Sqr(gp)*
      Sqr(Qq)*(Yu*mq2*Yu.adjoint()) + 20*Sqr(gp)*Sqr(Qu)*(Yu*mq2*Yu.adjoint())
      + 15*traceYuAdjYu*(Yu*Yu.adjoint()*mu2) + 5*AbsSqr(Lambdax)*(Yu*
      Yu.adjoint()*mu2) + Sqr(g1)*(Yu*Yu.adjoint()*mu2) - 15*Sqr(g2)*(Yu*
      Yu.adjoint()*mu2) - 10*Sqr(gp)*Sqr(QHu)*(Yu*Yu.adjoint()*mu2) - 10*Sqr(gp
      )*Sqr(Qq)*(Yu*Yu.adjoint()*mu2) + 10*Sqr(gp)*Sqr(Qu)*(Yu*Yu.adjoint()*mu2
      ) + 10*mHd2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()) + 10*mHu2*(Yu*Yd.adjoint()*
      Yd*Yu.adjoint()) + 10*(Yu*Yd.adjoint()*TYd*(TYu).adjoint()) + 20*mHu2*(Yu
      *Yu.adjoint()*Yu*Yu.adjoint()) + 10*(Yu*Yu.adjoint()*TYu*(TYu).adjoint())
      + 10*(Yu*(TYd).adjoint()*TYd*Yu.adjoint()) + 10*(Yu*(TYu).adjoint()*TYu*
      Yu.adjoint()) + 10*(TYu*Yd.adjoint()*Yd*(TYu).adjoint()) + 10*(TYu*
      Yu.adjoint()*Yu*(TYu).adjoint()) + 10*(TYu*(TYd).adjoint()*Yd*Yu.adjoint(
      )) + 10*(TYu*(TYu).adjoint()*Yu*Yu.adjoint()) + 5*(mu2*Yu*Yd.adjoint()*Yd
      *Yu.adjoint()) + 5*(mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint()) + 10*(Yu*mq2*
      Yd.adjoint()*Yd*Yu.adjoint()) + 10*(Yu*mq2*Yu.adjoint()*Yu*Yu.adjoint())
      + 10*(Yu*Yd.adjoint()*md2*Yd*Yu.adjoint()) + 10*(Yu*Yd.adjoint()*Yd*mq2*
      Yu.adjoint()) + 5*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*mu2) + 10*(Yu*
      Yu.adjoint()*mu2*Yu*Yu.adjoint()) + 10*(Yu*Yu.adjoint()*Yu*mq2*Yu.adjoint
      ()) + 5*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*mu2)) - 4*(20*Power(g3,4)*Tr23 -
      7.745966692414834*g1*(gp*Qu*(Tr2U114 + Tr2U141) + Tr31) + 15*gp*Qu*(gp*
      Qu*Tr2U144 + Tr34) + 4*Tr2U111*Sqr(g1))*UNITMATRIX(3) - 4*Conj(MassU)*Sqr
      (gp)*(15*(Sqr(QHu) + Sqr(Qq) - Sqr(Qu))*(2*MassU*(Yu*Yu.adjoint()) - TYu*
      Yu.adjoint()) + Qu*(-2*(MassB + 2*MassU)*(9*Qd + 9*Qe - 3*QHd + 3*QHu - 9
      *Ql + 9*Qq - 22*Qu)*Sqr(g1) + 5*Qu*(8*(MassG + 2*MassU)*Sqr(g3) + 9*MassU
      *Sqr(gp)*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) +
      18*Sqr(Qq) + Sqr(Qs) + 11*Sqr(Qu))))*UNITMATRIX(3)))))).real();


   return beta_mu2;
}

/**
 * Calculates the three-loop beta function of mu2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return beta_mu2;
}

} // namespace flexiblesusy
