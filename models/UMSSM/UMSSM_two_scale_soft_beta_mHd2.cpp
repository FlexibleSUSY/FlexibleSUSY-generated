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

// File generated at Sun 18 Oct 2015 12:18:49

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mHd2.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHd2;

   beta_mHd2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 + 2*gp*QHd*
      Tr14 + 6*traceconjTYdTpTYd + 2*traceconjTYeTpTYe + 6*tracemd2YdAdjYd + 2*
      traceme2YeAdjYe + 2*traceml2AdjYeYe + 6*tracemq2AdjYdYd + 6*mHd2*
      traceYdAdjYd + 2*mHd2*traceYeAdjYe + 2*mHd2*AbsSqr(Lambdax) + 2*mHu2*
      AbsSqr(Lambdax) + 2*ms2*AbsSqr(Lambdax) + 2*AbsSqr(TLambdax) - 1.2*AbsSqr
      (MassB)*Sqr(g1) - 6*AbsSqr(MassWB)*Sqr(g2) - 8*AbsSqr(MassU)*Sqr(gp)*Sqr(
      QHd)));


   return beta_mHd2;
}

/**
 * Calculates the two-loop beta function of mHd2.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qs = INPUT(Qs);
   const auto Qd = INPUT(Qd);
   const auto Qq = INPUT(Qq);
   const auto Qe = INPUT(Qe);
   const auto Ql = INPUT(Ql);
   const auto Qu = INPUT(Qu);
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd =
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd =
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe =
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe =
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double tracemd2YdAdjYdYdAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe =
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe =
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHd2;

   beta_mHd2 = Re(twoLoop*(6*Power(g2,4)*Tr22 - 3.0983866769659336*g1*gp*
      QHd*Tr2U114 - 3.0983866769659336*g1*gp*QHd*Tr2U141 - 3.0983866769659336*
      g1*Tr31 + 8*gp*QHd*Tr34 - 36*tracemd2YdAdjYdYdAdjYd - 6*
      tracemd2YdAdjYuYuAdjYd - 12*traceme2YeAdjYeYeAdjYe - 12*
      traceml2AdjYeYeAdjYeYe - 36*tracemq2AdjYdYdAdjYdYd - 6*
      tracemq2AdjYdYdAdjYuYu - 6*tracemq2AdjYuYuAdjYdYd - 6*
      tracemu2YuAdjYdYdAdjYu - 36*traceYdAdjTYdTYdAdjYd - 6*
      traceYdAdjTYuTYuAdjYd - 36*traceYdAdjYdTYdAdjTYd - 36*mHd2*
      traceYdAdjYdYdAdjYd - 6*traceYdAdjYuTYuAdjTYd - 6*mHd2*
      traceYdAdjYuYuAdjYd - 6*mHu2*traceYdAdjYuYuAdjYd - 12*
      traceYeAdjTYeTYeAdjYe - 12*traceYeAdjYeTYeAdjTYe - 12*mHd2*
      traceYeAdjYeYeAdjYe - 6*traceYuAdjTYdTYdAdjYu - 6*traceYuAdjYdTYdAdjTYu +
      33*Power(g2,4)*AbsSqr(MassWB) - 6*traceconjTYuTpTYu*AbsSqr(Lambdax) - 6*
      tracemq2AdjYuYu*AbsSqr(Lambdax) - 6*tracemu2YuAdjYu*AbsSqr(Lambdax) - 6*
      mHd2*traceYuAdjYu*AbsSqr(Lambdax) - 12*mHu2*traceYuAdjYu*AbsSqr(Lambdax)
      - 6*ms2*traceYuAdjYu*AbsSqr(Lambdax) - 6*traceYuAdjYu*AbsSqr(TLambdax) -
      24*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 6*traceAdjYuTYu*Conj(TLambdax)*
      Lambdax + 1.2*Tr2U111*Sqr(g1) - 0.8*traceconjTYdTpTYd*Sqr(g1) + 0.8*MassB
      *traceconjTYdTpYd*Sqr(g1) + 2.4*traceconjTYeTpTYe*Sqr(g1) - 2.4*MassB*
      traceconjTYeTpYe*Sqr(g1) - 0.8*tracemd2YdAdjYd*Sqr(g1) + 2.4*
      traceme2YeAdjYe*Sqr(g1) + 2.4*traceml2AdjYeYe*Sqr(g1) - 0.8*
      tracemq2AdjYdYd*Sqr(g1) - 0.8*mHd2*traceYdAdjYd*Sqr(g1) + 2.4*mHd2*
      traceYeAdjYe*Sqr(g1) + 3.6*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 1.8*MassB*
      Conj(MassWB)*Sqr(g1)*Sqr(g2) + 32*traceconjTYdTpTYd*Sqr(g3) - 32*MassG*
      traceconjTYdTpYd*Sqr(g3) + 32*tracemd2YdAdjYd*Sqr(g3) + 32*
      tracemq2AdjYdYd*Sqr(g3) + 32*mHd2*traceYdAdjYd*Sqr(g3) + 64*traceYdAdjYd*
      AbsSqr(MassG)*Sqr(g3) - 32*traceAdjYdTYd*Conj(MassG)*Sqr(g3) + 12*
      traceconjTYdTpTYd*Sqr(gp)*Sqr(Qd) - 12*MassU*traceconjTYdTpYd*Sqr(gp)*Sqr
      (Qd) + 12*tracemd2YdAdjYd*Sqr(gp)*Sqr(Qd) + 12*tracemq2AdjYdYd*Sqr(gp)*
      Sqr(Qd) + 12*mHd2*traceYdAdjYd*Sqr(gp)*Sqr(Qd) + 4*traceconjTYeTpTYe*Sqr(
      gp)*Sqr(Qe) - 4*MassU*traceconjTYeTpYe*Sqr(gp)*Sqr(Qe) + 4*
      traceme2YeAdjYe*Sqr(gp)*Sqr(Qe) + 4*traceml2AdjYeYe*Sqr(gp)*Sqr(Qe) + 4*
      mHd2*traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 8*Tr2U144*Sqr(gp)*Sqr(QHd) - 12*
      traceconjTYdTpTYd*Sqr(gp)*Sqr(QHd) + 12*MassU*traceconjTYdTpYd*Sqr(gp)*
      Sqr(QHd) - 4*traceconjTYeTpTYe*Sqr(gp)*Sqr(QHd) + 4*MassU*
      traceconjTYeTpYe*Sqr(gp)*Sqr(QHd) - 12*tracemd2YdAdjYd*Sqr(gp)*Sqr(QHd) -
      4*traceme2YeAdjYe*Sqr(gp)*Sqr(QHd) - 4*traceml2AdjYeYe*Sqr(gp)*Sqr(QHd)
      - 12*tracemq2AdjYdYd*Sqr(gp)*Sqr(QHd) - 12*mHd2*traceYdAdjYd*Sqr(gp)*Sqr(
      QHd) - 4*mHd2*traceYeAdjYe*Sqr(gp)*Sqr(QHd) - 4*mHd2*AbsSqr(Lambdax)*Sqr(
      gp)*Sqr(QHd) - 4*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) - 4*ms2*AbsSqr(
      Lambdax)*Sqr(gp)*Sqr(QHd) - 4*AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHd) + 4*MassU
      *Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(QHd) + 24*AbsSqr(MassWB)*Sqr(g2)*Sqr(
      gp)*Sqr(QHd) + 12*MassU*Conj(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 0.04*Conj
      (MassB)*Sqr(g1)*(20*traceAdjYdTYd - 60*traceAdjYeTYe - 40*MassB*
      traceYdAdjYd + 120*MassB*traceYeAdjYe + 621*MassB*Sqr(g1) + 90*MassB*Sqr(
      g2) + 45*MassWB*Sqr(g2) - 360*MassB*Qd*QHd*Sqr(gp) - 180*MassU*Qd*QHd*Sqr
      (gp) - 360*MassB*Qe*QHd*Sqr(gp) - 180*MassU*Qe*QHd*Sqr(gp) - 120*MassB*
      QHd*QHu*Sqr(gp) - 60*MassU*QHd*QHu*Sqr(gp) + 360*MassB*QHd*Ql*Sqr(gp) +
      180*MassU*QHd*Ql*Sqr(gp) - 360*MassB*QHd*Qq*Sqr(gp) - 180*MassU*QHd*Qq*
      Sqr(gp) + 720*MassB*QHd*Qu*Sqr(gp) + 360*MassU*QHd*Qu*Sqr(gp) + 240*MassB
      *Sqr(gp)*Sqr(QHd) + 120*MassU*Sqr(gp)*Sqr(QHd)) + 4*mHd2*AbsSqr(Lambdax)*
      Sqr(gp)*Sqr(QHu) + 4*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 4*ms2*AbsSqr
      (Lambdax)*Sqr(gp)*Sqr(QHu) + 4*AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHu) - 4*
      MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(QHu) + 4*traceconjTYeTpTYe*Sqr(
      gp)*Sqr(Ql) - 4*MassU*traceconjTYeTpYe*Sqr(gp)*Sqr(Ql) + 4*
      traceme2YeAdjYe*Sqr(gp)*Sqr(Ql) + 4*traceml2AdjYeYe*Sqr(gp)*Sqr(Ql) + 4*
      mHd2*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 12*traceconjTYdTpTYd*Sqr(gp)*Sqr(Qq)
      - 12*MassU*traceconjTYdTpYd*Sqr(gp)*Sqr(Qq) + 12*tracemd2YdAdjYd*Sqr(gp)*
      Sqr(Qq) + 12*tracemq2AdjYdYd*Sqr(gp)*Sqr(Qq) + 12*mHd2*traceYdAdjYd*Sqr(
      gp)*Sqr(Qq) + 4*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 4*mHu2*AbsSqr(
      Lambdax)*Sqr(gp)*Sqr(Qs) + 4*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 4*
      AbsSqr(TLambdax)*Sqr(gp)*Sqr(Qs) - 4*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)
      *Sqr(Qs) - 12*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 12*mHu2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) - 12*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 0.8*
      Conj(MassU)*Sqr(gp)*(-9*MassB*Qd*QHd*Sqr(g1) - 18*MassU*Qd*QHd*Sqr(g1) -
      9*MassB*Qe*QHd*Sqr(g1) - 18*MassU*Qe*QHd*Sqr(g1) - 3*MassB*QHd*QHu*Sqr(g1
      ) - 6*MassU*QHd*QHu*Sqr(g1) + 9*MassB*QHd*Ql*Sqr(g1) + 18*MassU*QHd*Ql*
      Sqr(g1) - 9*MassB*QHd*Qq*Sqr(g1) - 18*MassU*QHd*Qq*Sqr(g1) + 18*MassB*QHd
      *Qu*Sqr(g1) + 36*MassU*QHd*Qu*Sqr(g1) + 120*MassU*Power(QHd,4)*Sqr(gp) -
      15*traceAdjYdTYd*Sqr(Qd) - 5*traceAdjYeTYe*Sqr(Qe) + 10*MassU*
      traceYeAdjYe*Sqr(Qe) + 15*traceAdjYdTYd*Sqr(QHd) + 5*traceAdjYeTYe*Sqr(
      QHd) - 10*MassU*traceYeAdjYe*Sqr(QHd) + 6*MassB*Sqr(g1)*Sqr(QHd) + 12*
      MassU*Sqr(g1)*Sqr(QHd) + 30*MassU*Sqr(g2)*Sqr(QHd) + 15*MassWB*Sqr(g2)*
      Sqr(QHd) + 270*MassU*Sqr(gp)*Sqr(Qd)*Sqr(QHd) + 90*MassU*Sqr(gp)*Sqr(Qe)*
      Sqr(QHd) + 60*MassU*Sqr(gp)*Sqr(QHd)*Sqr(QHu) - 5*traceAdjYeTYe*Sqr(Ql) +
      10*MassU*traceYeAdjYe*Sqr(Ql) + 180*MassU*Sqr(gp)*Sqr(QHd)*Sqr(Ql) - 15*
      traceAdjYdTYd*Sqr(Qq) + 540*MassU*Sqr(gp)*Sqr(QHd)*Sqr(Qq) + 30*MassU*
      traceYdAdjYd*(Sqr(Qd) - Sqr(QHd) + Sqr(Qq)) + 30*MassU*Sqr(gp)*Sqr(QHd)*
      Sqr(Qs) + 270*MassU*Sqr(gp)*Sqr(QHd)*Sqr(Qu) - 5*Conj(Lambdax)*(Sqr(QHd)
      - Sqr(QHu) - Sqr(Qs))*(2*MassU*Lambdax - TLambdax)) - 6*traceconjTYuTpYu*
      Conj(Lambdax)*TLambdax));


   return beta_mHd2;
}

/**
 * Calculates the three-loop beta function of mHd2.
 *
 * @return three-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

} // namespace flexiblesusy
