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

// File generated at Mon 23 Feb 2015 13:04:21

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHu2.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHu = INPUT(QHu);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHu2;

   beta_mHu2 = oneOver16PiSqr*(0.7745966692414834*g1*Tr11 + 2*gp*QHu*Tr14
      + 6*traceconjTYuTpTYu + 6*tracemq2AdjYuYu + 6*tracemu2YuAdjYu + 6*mHu2*
      traceYuAdjYu + 2*mHd2*AbsSqr(Lambdax) + 2*mHu2*AbsSqr(Lambdax) + 2*ms2*
      AbsSqr(Lambdax) + 2*AbsSqr(TLambdax) - 1.2*AbsSqr(MassB)*Sqr(g1) - 6*
      AbsSqr(MassWB)*Sqr(g2) - 8*AbsSqr(MassU)*Sqr(gp)*Sqr(QHu));


   return beta_mHu2;
}

/**
 * Calculates the two-loop beta function of mHu2.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto QHd = INPUT(QHd);
   const auto Qs = INPUT(Qs);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto Ql = INPUT(Ql);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu =
      TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu =
      TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHu2;

   beta_mHu2 = twoLoop*(6*Power(g2,4)*Tr22 + 3.0983866769659336*g1*gp*QHu
      *Tr2U114 + 3.0983866769659336*g1*gp*QHu*Tr2U141 + 3.0983866769659336*g1*
      Tr31 + 8*gp*QHu*Tr34 - 6*tracemd2YdAdjYuYuAdjYd - 6*
      tracemq2AdjYdYdAdjYuYu - 6*tracemq2AdjYuYuAdjYdYd - 36*
      tracemq2AdjYuYuAdjYuYu - 6*tracemu2YuAdjYdYdAdjYu - 36*
      tracemu2YuAdjYuYuAdjYu - 6*traceYdAdjTYuTYuAdjYd - 6*
      traceYdAdjYuTYuAdjTYd - 6*mHd2*traceYdAdjYuYuAdjYd - 6*mHu2*
      traceYdAdjYuYuAdjYd - 6*traceYuAdjTYdTYdAdjYu - 36*traceYuAdjTYuTYuAdjYu
      - 6*traceYuAdjYdTYdAdjTYu - 36*traceYuAdjYuTYuAdjTYu - 36*mHu2*
      traceYuAdjYuYuAdjYu + 33*Power(g2,4)*AbsSqr(MassWB) - 6*traceconjTYdTpTYd
      *AbsSqr(Lambdax) - 2*traceconjTYeTpTYe*AbsSqr(Lambdax) - 6*
      tracemd2YdAdjYd*AbsSqr(Lambdax) - 2*traceme2YeAdjYe*AbsSqr(Lambdax) - 2*
      traceml2AdjYeYe*AbsSqr(Lambdax) - 6*tracemq2AdjYdYd*AbsSqr(Lambdax) - 12*
      mHd2*traceYdAdjYd*AbsSqr(Lambdax) - 6*mHu2*traceYdAdjYd*AbsSqr(Lambdax) -
      6*ms2*traceYdAdjYd*AbsSqr(Lambdax) - 4*mHd2*traceYeAdjYe*AbsSqr(Lambdax)
      - 2*mHu2*traceYeAdjYe*AbsSqr(Lambdax) - 2*ms2*traceYeAdjYe*AbsSqr(
      Lambdax) - 6*traceYdAdjYd*AbsSqr(TLambdax) - 2*traceYeAdjYe*AbsSqr(
      TLambdax) - 24*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 6*traceAdjYdTYd*Conj(
      TLambdax)*Lambdax - 2*traceAdjYeTYe*Conj(TLambdax)*Lambdax + 1.2*Tr2U111*
      Sqr(g1) + 1.6*traceconjTYuTpTYu*Sqr(g1) - 1.6*MassB*traceconjTYuTpYu*Sqr(
      g1) + 1.6*tracemq2AdjYuYu*Sqr(g1) + 1.6*tracemu2YuAdjYu*Sqr(g1) + 1.6*
      mHu2*traceYuAdjYu*Sqr(g1) + 3.6*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 1.8*
      MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 32*traceconjTYuTpTYu*Sqr(g3) - 32*
      MassG*traceconjTYuTpYu*Sqr(g3) + 32*tracemq2AdjYuYu*Sqr(g3) + 32*
      tracemu2YuAdjYu*Sqr(g3) + 32*mHu2*traceYuAdjYu*Sqr(g3) + 64*traceYuAdjYu*
      AbsSqr(MassG)*Sqr(g3) - 32*traceAdjYuTYu*Conj(MassG)*Sqr(g3) + 0.04*Conj(
      MassB)*Sqr(g1)*(-40*traceAdjYuTYu + 80*MassB*traceYuAdjYu + 3*(207*MassB*
      Sqr(g1) + 5*(3*(2*MassB + MassWB)*Sqr(g2) + 4*(2*MassB + MassU)*QHu*(3*Qd
      + 3*Qe - QHd + 2*QHu - 3*Ql + 3*Qq - 6*Qu)*Sqr(gp)))) + 4*mHd2*AbsSqr(
      Lambdax)*Sqr(gp)*Sqr(QHd) + 4*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 4*
      ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 4*AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHd
      ) - 4*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(QHd) + 8*Tr2U144*Sqr(gp)*
      Sqr(QHu) - 12*traceconjTYuTpTYu*Sqr(gp)*Sqr(QHu) + 12*MassU*
      traceconjTYuTpYu*Sqr(gp)*Sqr(QHu) - 12*tracemq2AdjYuYu*Sqr(gp)*Sqr(QHu) -
      12*tracemu2YuAdjYu*Sqr(gp)*Sqr(QHu) - 12*mHu2*traceYuAdjYu*Sqr(gp)*Sqr(
      QHu) - 4*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) - 4*mHu2*AbsSqr(Lambdax)*
      Sqr(gp)*Sqr(QHu) - 4*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) - 4*AbsSqr(
      TLambdax)*Sqr(gp)*Sqr(QHu) + 4*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(
      QHu) + 24*AbsSqr(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 12*MassU*Conj(MassWB)
      *Sqr(g2)*Sqr(gp)*Sqr(QHu) + 12*traceconjTYuTpTYu*Sqr(gp)*Sqr(Qq) - 12*
      MassU*traceconjTYuTpYu*Sqr(gp)*Sqr(Qq) + 12*tracemq2AdjYuYu*Sqr(gp)*Sqr(
      Qq) + 12*tracemu2YuAdjYu*Sqr(gp)*Sqr(Qq) + 12*mHu2*traceYuAdjYu*Sqr(gp)*
      Sqr(Qq) + 4*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 4*mHu2*AbsSqr(Lambdax)
      *Sqr(gp)*Sqr(Qs) + 4*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 4*AbsSqr(
      TLambdax)*Sqr(gp)*Sqr(Qs) - 4*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(Qs
      ) + 12*traceconjTYuTpTYu*Sqr(gp)*Sqr(Qu) - 12*MassU*traceconjTYuTpYu*Sqr(
      gp)*Sqr(Qu) + 12*tracemq2AdjYuYu*Sqr(gp)*Sqr(Qu) + 12*tracemu2YuAdjYu*Sqr
      (gp)*Sqr(Qu) + 12*mHu2*traceYuAdjYu*Sqr(gp)*Sqr(Qu) - 12*mHd2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) - 12*mHu2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 12*ms2
      *Sqr(Conj(Lambdax))*Sqr(Lambdax) + 0.8*Conj(MassU)*Sqr(gp)*(3*(5*
      traceAdjYuTYu*(Sqr(QHu) - Sqr(Qq) - Sqr(Qu)) - 10*MassU*traceYuAdjYu*(Sqr
      (QHu) - Sqr(Qq) - Sqr(Qu)) + QHu*((MassB + 2*MassU)*(3*Qd + 3*Qe - QHd +
      2*QHu - 3*Ql + 3*Qq - 6*Qu)*Sqr(g1) + 5*QHu*((2*MassU + MassWB)*Sqr(g2) +
      2*MassU*Sqr(gp)*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(QHd) + 4*Sqr(QHu) + 6*Sqr
      (Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu))))) + 5*Conj(Lambdax)*(Sqr(QHd) -
      Sqr(QHu) + Sqr(Qs))*(2*MassU*Lambdax - TLambdax)) - 6*traceconjTYdTpYd*
      Conj(Lambdax)*TLambdax - 2*traceconjTYeTpYe*Conj(Lambdax)*TLambdax);


   return beta_mHu2;
}

} // namespace flexiblesusy
