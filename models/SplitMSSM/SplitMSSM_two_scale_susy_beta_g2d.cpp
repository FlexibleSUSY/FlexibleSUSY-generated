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

// File generated at Mon 19 Sep 2016 09:48:20

#include "SplitMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g2d.
 *
 * @return one-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2d;

   beta_g2d = Re(oneOver16PiSqr*(2.75*Power(g2d,3) + g2u*gYd*gYu + 3*g2d*
      traceYdAdjYd + g2d*traceYeAdjYe + 3*g2d*traceYuAdjYu - 0.45*g2d*Sqr(g1) -
      8.25*g2d*Sqr(g2) + g2d*Sqr(g2u) + 0.75*g2d*Sqr(gYd) + 0.5*g2d*Sqr(gYu)))
      ;


   return beta_g2d;
}

/**
 * Calculates the two-loop beta function of g2d.
 *
 * @return two-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2d;

   beta_g2d = Re(twoLoop*(0.585*Power(g1,4)*g2d - 34.083333333333336*
      Power(g2,4)*g2d - 3.5*Power(g2d,5) - 1.375*g2d*Power(g2u,4) - 0.3125*g2d*
      Power(gYd,4) - 2.25*Power(g2u,3)*gYd*gYu - 1.5*g2u*Power(gYd,3)*gYu -
      1.25*g2u*gYd*Power(gYu,3) - 0.5625*g2d*Power(gYu,4) - 6.75*g2d*
      traceYdAdjYdYdAdjYd + 1.5*g2d*traceYdAdjYuYuAdjYd - 2.25*g2d*
      traceYeAdjYeYeAdjYe - 5.625*Power(g2d,3)*traceYuAdjYu - 3*g2u*gYd*gYu*
      traceYuAdjYu - 6.75*g2d*traceYuAdjYuYuAdjYu - 5*Power(g2d,3)*Lambdax -
      g2u*gYd*gYu*Lambdax + 2.71875*Power(g2d,3)*Sqr(g1) + 0.15*g2u*gYd*gYu*Sqr
      (g1) + 2.125*g2d*traceYuAdjYu*Sqr(g1) + 27.34375*Power(g2d,3)*Sqr(g2) +
      2.25*g2u*gYd*gYu*Sqr(g2) + 5.625*g2d*traceYuAdjYu*Sqr(g2) + 0.45*g2d*Sqr(
      g1)*Sqr(g2) - 4*g2u*gYd*gYu*Sqr(g2d) - 3.375*Power(g2d,3)*Sqr(g2u) + 0.75
      *g2d*traceYuAdjYu*Sqr(g2u) - g2d*Lambdax*Sqr(g2u) + 0.15*g2d*Sqr(g1)*Sqr(
      g2u) + 4.25*g2d*Sqr(g2)*Sqr(g2u) + 20*g2d*traceYuAdjYu*Sqr(g3) - 3.6875*
      Power(g2d,3)*Sqr(gYd) - 1.125*g2d*traceYuAdjYu*Sqr(gYd) - g2d*Lambdax*Sqr
      (gYd) + 0.39375*g2d*Sqr(g1)*Sqr(gYd) + 3.46875*g2d*Sqr(g2)*Sqr(gYd) -
      1.9375*g2d*Sqr(g2u)*Sqr(gYd) + 0.125*traceYdAdjYd*(-45*Power(g2d,3) - 24*
      g2u*gYd*gYu + 5*g2d*Sqr(g1) + 45*g2d*Sqr(g2) + 6*g2d*Sqr(g2u) + 160*g2d*
      Sqr(g3) - 9*g2d*Sqr(gYd)) + 0.125*traceYeAdjYe*(-15*Power(g2d,3) - 8*g2u*
      gYd*gYu + 15*g2d*Sqr(g1) + 15*g2d*Sqr(g2) + 2*g2d*Sqr(g2u) - 3*g2d*Sqr(
      gYd)) - 0.9375*Power(g2d,3)*Sqr(gYu) + 0.1875*g2d*Sqr(g1)*Sqr(gYu) +
      0.9375*g2d*Sqr(g2)*Sqr(gYu) - 0.8125*g2d*Sqr(g2u)*Sqr(gYu) - 1.5*g2d*Sqr(
      gYd)*Sqr(gYu) + 1.5*g2d*Sqr(Lambdax)));


   return beta_g2d;
}

/**
 * Calculates the three-loop beta function of g2d.
 *
 * @return three-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2d;

   beta_g2d = 0;


   return beta_g2d;
}

} // namespace flexiblesusy
