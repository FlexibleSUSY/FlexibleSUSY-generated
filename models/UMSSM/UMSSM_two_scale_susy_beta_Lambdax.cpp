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

// File generated at Tue 7 Jul 2015 12:46:26

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambdax.
 *
 * @return one-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_Lambdax_one_loop(const Susy_traces& susy_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qs = INPUT(Qs);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(oneOver16PiSqr*(3*traceYdAdjYd*Lambdax +
      traceYeAdjYe*Lambdax + 3*traceYuAdjYu*Lambdax - 0.6*Lambdax*Sqr(g1) - 3*
      Lambdax*Sqr(g2) - 2*Lambdax*Sqr(gp)*Sqr(QHd) - 2*Lambdax*Sqr(gp)*Sqr(QHu)
      - 2*Lambdax*Sqr(gp)*Sqr(Qs) + 4*Conj(Lambdax)*Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the two-loop beta function of Lambdax.
 *
 * @return two-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_Lambdax_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qe = INPUT(Qe);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(-0.02*twoLoop*Lambdax*(-207*Power(g1,4) - 375*Power(
      g2,4) - 400*Power(gp,4)*Power(QHd,4) - 400*Power(gp,4)*Power(QHu,4) - 300
      *Power(gp,4)*Power(Qs,4) + 450*traceYdAdjYdYdAdjYd + 300*
      traceYdAdjYuYuAdjYd + 150*traceYeAdjYeYeAdjYe + 450*traceYuAdjYuYuAdjYu -
      60*traceYeAdjYe*Sqr(g1) - 40*traceYuAdjYu*Sqr(g1) - 90*Sqr(g1)*Sqr(g2) -
      800*traceYuAdjYu*Sqr(g3) + 180*Qd*QHd*Sqr(g1)*Sqr(gp) + 180*Qe*QHd*Sqr(
      g1)*Sqr(gp) - 180*Qd*QHu*Sqr(g1)*Sqr(gp) - 180*Qe*QHu*Sqr(g1)*Sqr(gp) +
      120*QHd*QHu*Sqr(g1)*Sqr(gp) - 180*QHd*Ql*Sqr(g1)*Sqr(gp) + 180*QHu*Ql*Sqr
      (g1)*Sqr(gp) + 180*QHd*Qq*Sqr(g1)*Sqr(gp) - 180*QHu*Qq*Sqr(g1)*Sqr(gp) -
      360*QHd*Qu*Sqr(g1)*Sqr(gp) + 360*QHu*Qu*Sqr(g1)*Sqr(gp) - 100*
      traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 100*traceYeAdjYe*Sqr(gp)*Sqr(QHd) - 120*
      Sqr(g1)*Sqr(gp)*Sqr(QHd) - 300*Sqr(g2)*Sqr(gp)*Sqr(QHd) - 900*Power(gp,4)
      *Sqr(Qd)*Sqr(QHd) - 300*Power(gp,4)*Sqr(Qe)*Sqr(QHd) + 300*traceYuAdjYu*
      Sqr(gp)*Sqr(QHu) - 120*Sqr(g1)*Sqr(gp)*Sqr(QHu) - 300*Sqr(g2)*Sqr(gp)*Sqr
      (QHu) - 900*Power(gp,4)*Sqr(Qd)*Sqr(QHu) - 300*Power(gp,4)*Sqr(Qe)*Sqr(
      QHu) - 400*Power(gp,4)*Sqr(QHd)*Sqr(QHu) - 10*AbsSqr(Lambdax)*(-45*
      traceYdAdjYd - 15*traceYeAdjYe - 45*traceYuAdjYu + 6*Sqr(g1) + 30*Sqr(g2)
      + 20*Sqr(gp)*Sqr(QHd) + 20*Sqr(gp)*Sqr(QHu)) - 100*traceYeAdjYe*Sqr(gp)*
      Sqr(Ql) - 600*Power(gp,4)*Sqr(QHd)*Sqr(Ql) - 600*Power(gp,4)*Sqr(QHu)*Sqr
      (Ql) - 300*traceYuAdjYu*Sqr(gp)*Sqr(Qq) - 1800*Power(gp,4)*Sqr(QHd)*Sqr(
      Qq) - 1800*Power(gp,4)*Sqr(QHu)*Sqr(Qq) + 20*traceYdAdjYd*(Sqr(g1) - 5*(8
      *Sqr(g3) + 3*Sqr(gp)*(Sqr(Qd) - Sqr(QHd) + Sqr(Qq)))) - 900*Power(gp,4)*
      Sqr(Qd)*Sqr(Qs) - 300*Power(gp,4)*Sqr(Qe)*Sqr(Qs) - 300*Power(gp,4)*Sqr(
      QHd)*Sqr(Qs) - 300*Power(gp,4)*Sqr(QHu)*Sqr(Qs) - 600*Power(gp,4)*Sqr(Ql)
      *Sqr(Qs) - 1800*Power(gp,4)*Sqr(Qq)*Sqr(Qs) - 300*traceYuAdjYu*Sqr(gp)*
      Sqr(Qu) - 900*Power(gp,4)*Sqr(QHd)*Sqr(Qu) - 900*Power(gp,4)*Sqr(QHu)*Sqr
      (Qu) - 900*Power(gp,4)*Sqr(Qs)*Sqr(Qu) + 500*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the three-loop beta function of Lambdax.
 *
 * @return three-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_Lambdax_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
