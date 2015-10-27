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

// File generated at Tue 27 Oct 2015 15:07:11

#include "SplitMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of gYd.
 *
 * @return one-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_gYd;

   beta_gYd = Re(oneOver16PiSqr*(1.25*Power(gYd,3) + 3*g2d*g2u*gYu + 3*
      gYd*traceYdAdjYd + gYd*traceYeAdjYe + 3*gYd*traceYuAdjYu - 0.45*gYd*Sqr(
      g1) - 2.25*gYd*Sqr(g2) + 2.25*gYd*Sqr(g2d) + 1.5*gYd*Sqr(g2u) + 2*gYd*Sqr
      (gYu)));


   return beta_gYd;
}

/**
 * Calculates the two-loop beta function of gYd.
 *
 * @return two-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_gYd;

   beta_gYd = Re(twoLoop*(0.585*Power(g1,4)*gYd - 4.25*Power(g2,4)*gYd -
      6.1875*Power(g2d,4)*gYd - 2.8125*Power(g2u,4)*gYd - 0.75*Power(gYd,5) -
      4.5*Power(g2d,3)*g2u*gYu - 8.25*g2d*Power(g2u,3)*gYu - 2.25*g2d*g2u*Power
      (gYu,3) - 2.25*gYd*Power(gYu,4) - 6.75*gYd*traceYdAdjYdYdAdjYd + 1.5*gYd*
      traceYdAdjYuYuAdjYd - 2.25*gYd*traceYeAdjYeYeAdjYe - 3.375*Power(gYd,3)*
      traceYuAdjYu - 9*g2d*g2u*gYu*traceYuAdjYu - 6.75*gYd*traceYuAdjYuYuAdjYu
      - 3*Power(gYd,3)*Lambdax - 3*g2d*g2u*gYu*Lambdax + 1.93125*Power(gYd,3)*
      Sqr(g1) + 0.45*g2d*g2u*gYu*Sqr(g1) + 2.125*gYd*traceYuAdjYu*Sqr(g1) +
      5.15625*Power(gYd,3)*Sqr(g2) + 12.75*g2d*g2u*gYu*Sqr(g2) + 5.625*gYd*
      traceYuAdjYu*Sqr(g2) - 1.35*gYd*Sqr(g1)*Sqr(g2) - 0.5625*Power(gYd,3)*Sqr
      (g2d) - 3.375*gYd*traceYuAdjYu*Sqr(g2d) - 3*gYd*Lambdax*Sqr(g2d) +
      1.18125*gYd*Sqr(g1)*Sqr(g2d) + 17.15625*gYd*Sqr(g2)*Sqr(g2d) - 1.6875*
      Power(gYd,3)*Sqr(g2u) + 0.5625*gYd*Sqr(g1)*Sqr(g2u) + 10.3125*gYd*Sqr(g2)
      *Sqr(g2u) - 2.625*gYd*Sqr(g2d)*Sqr(g2u) + 20*gYd*traceYuAdjYu*Sqr(g3) - 6
      *g2d*g2u*gYu*Sqr(gYd) - 3.75*Power(gYd,3)*Sqr(gYu) - 5.25*gYd*
      traceYuAdjYu*Sqr(gYu) - 3*gYd*Lambdax*Sqr(gYu) + 0.075*gYd*Sqr(g1)*Sqr(
      gYu) + 4.875*gYd*Sqr(g2)*Sqr(gYu) - 4.6875*gYd*Sqr(g2d)*Sqr(gYu) - 4.6875
      *gYd*Sqr(g2u)*Sqr(gYu) + 0.125*traceYdAdjYd*(-27*Power(gYd,3) - 72*g2d*
      g2u*gYu + 5*gYd*Sqr(g1) + 45*gYd*Sqr(g2) - 27*gYd*Sqr(g2d) + 160*gYd*Sqr(
      g3) - 42*gYd*Sqr(gYu)) + 0.125*traceYeAdjYe*(-9*Power(gYd,3) - 24*g2d*g2u
      *gYu + 15*gYd*Sqr(g1) + 15*gYd*Sqr(g2) - 9*gYd*Sqr(g2d) - 14*gYd*Sqr(gYu)
      ) + 1.5*gYd*Sqr(Lambdax)));


   return beta_gYd;
}

/**
 * Calculates the three-loop beta function of gYd.
 *
 * @return three-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYd;

   beta_gYd = 0;


   return beta_gYd;
}

} // namespace flexiblesusy
