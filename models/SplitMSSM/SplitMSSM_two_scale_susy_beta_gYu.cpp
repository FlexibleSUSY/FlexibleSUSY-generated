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

// File generated at Tue 27 Oct 2015 15:07:12

#include "SplitMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of gYu.
 *
 * @return one-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_gYu;

   beta_gYu = Re(oneOver16PiSqr*(3*g2d*g2u*gYd + 1.25*Power(gYu,3) + 3*
      gYu*traceYdAdjYd + gYu*traceYeAdjYe + 3*gYu*traceYuAdjYu - 0.45*gYu*Sqr(
      g1) - 2.25*gYu*Sqr(g2) + 1.5*gYu*Sqr(g2d) + 2.25*gYu*Sqr(g2u) + 2*gYu*Sqr
      (gYd)));


   return beta_gYu;
}

/**
 * Calculates the two-loop beta function of gYu.
 *
 * @return two-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_gYu;

   beta_gYu = Re(twoLoop*(-8.25*Power(g2d,3)*g2u*gYd - 4.5*g2d*Power(g2u,
      3)*gYd - 2.25*g2d*g2u*Power(gYd,3) + 0.585*Power(g1,4)*gYu - 4.25*Power(
      g2,4)*gYu - 2.8125*Power(g2d,4)*gYu - 6.1875*Power(g2u,4)*gYu - 2.25*
      Power(gYd,4)*gYu - 0.75*Power(gYu,5) - 6.75*gYu*traceYdAdjYdYdAdjYd + 1.5
      *gYu*traceYdAdjYuYuAdjYd - 2.25*gYu*traceYeAdjYeYeAdjYe - 9*g2d*g2u*gYd*
      traceYuAdjYu - 3.375*Power(gYu,3)*traceYuAdjYu - 6.75*gYu*
      traceYuAdjYuYuAdjYu - 3*g2d*g2u*gYd*Lambdax - 3*Power(gYu,3)*Lambdax +
      0.45*g2d*g2u*gYd*Sqr(g1) + 1.93125*Power(gYu,3)*Sqr(g1) + 2.125*gYu*
      traceYuAdjYu*Sqr(g1) + 12.75*g2d*g2u*gYd*Sqr(g2) + 5.15625*Power(gYu,3)*
      Sqr(g2) + 5.625*gYu*traceYuAdjYu*Sqr(g2) - 1.35*gYu*Sqr(g1)*Sqr(g2) -
      1.6875*Power(gYu,3)*Sqr(g2d) + 0.5625*gYu*Sqr(g1)*Sqr(g2d) + 10.3125*gYu*
      Sqr(g2)*Sqr(g2d) - 0.5625*Power(gYu,3)*Sqr(g2u) - 3.375*gYu*traceYuAdjYu*
      Sqr(g2u) - 3*gYu*Lambdax*Sqr(g2u) + 1.18125*gYu*Sqr(g1)*Sqr(g2u) +
      17.15625*gYu*Sqr(g2)*Sqr(g2u) - 2.625*gYu*Sqr(g2d)*Sqr(g2u) + 20*gYu*
      traceYuAdjYu*Sqr(g3) - 3.75*Power(gYu,3)*Sqr(gYd) - 5.25*gYu*traceYuAdjYu
      *Sqr(gYd) - 3*gYu*Lambdax*Sqr(gYd) + 0.075*gYu*Sqr(g1)*Sqr(gYd) + 4.875*
      gYu*Sqr(g2)*Sqr(gYd) - 4.6875*gYu*Sqr(g2d)*Sqr(gYd) - 4.6875*gYu*Sqr(g2u)
      *Sqr(gYd) - 6*g2d*g2u*gYd*Sqr(gYu) - 0.125*traceYeAdjYe*(24*g2d*g2u*gYd +
      gYu*(-15*Sqr(g1) - 15*Sqr(g2) + 9*Sqr(g2u) + 14*Sqr(gYd) + 9*Sqr(gYu)))
      - 0.125*traceYdAdjYd*(72*g2d*g2u*gYd + gYu*(-5*Sqr(g1) - 45*Sqr(g2) + 27*
      Sqr(g2u) - 160*Sqr(g3) + 42*Sqr(gYd) + 27*Sqr(gYu))) + 1.5*gYu*Sqr(
      Lambdax)));


   return beta_gYu;
}

/**
 * Calculates the three-loop beta function of gYu.
 *
 * @return three-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYu;

   beta_gYu = 0;


   return beta_gYu;
}

} // namespace flexiblesusy
