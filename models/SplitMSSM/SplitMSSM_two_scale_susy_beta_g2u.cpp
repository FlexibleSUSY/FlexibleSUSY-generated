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

// File generated at Wed 29 Jun 2016 11:28:19

#include "SplitMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g2u.
 *
 * @return one-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(oneOver16PiSqr*(2.75*Power(g2u,3) + g2d*gYd*gYu + 3*g2u*
      traceYdAdjYd + g2u*traceYeAdjYe + 3*g2u*traceYuAdjYu - 0.45*g2u*Sqr(g1) -
      8.25*g2u*Sqr(g2) + g2u*Sqr(g2d) + 0.5*g2u*Sqr(gYd) + 0.75*g2u*Sqr(gYu)))
      ;


   return beta_g2u;
}

/**
 * Calculates the two-loop beta function of g2u.
 *
 * @return two-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(twoLoop*(0.585*Power(g1,4)*g2u - 34.083333333333336*
      Power(g2,4)*g2u - 1.375*Power(g2d,4)*g2u - 3.5*Power(g2u,5) - 0.5625*g2u*
      Power(gYd,4) - 2.25*Power(g2d,3)*gYd*gYu - 1.25*g2d*Power(gYd,3)*gYu -
      1.5*g2d*gYd*Power(gYu,3) - 0.3125*g2u*Power(gYu,4) - 6.75*g2u*
      traceYdAdjYdYdAdjYd + 1.5*g2u*traceYdAdjYuYuAdjYd - 2.25*g2u*
      traceYeAdjYeYeAdjYe - 5.625*Power(g2u,3)*traceYuAdjYu - 3*g2d*gYd*gYu*
      traceYuAdjYu - 6.75*g2u*traceYuAdjYuYuAdjYu - 5*Power(g2u,3)*Lambdax -
      g2d*gYd*gYu*Lambdax + 2.71875*Power(g2u,3)*Sqr(g1) + 0.15*g2d*gYd*gYu*Sqr
      (g1) + 2.125*g2u*traceYuAdjYu*Sqr(g1) + 27.34375*Power(g2u,3)*Sqr(g2) +
      2.25*g2d*gYd*gYu*Sqr(g2) + 5.625*g2u*traceYuAdjYu*Sqr(g2) + 0.45*g2u*Sqr(
      g1)*Sqr(g2) - 3.375*Power(g2u,3)*Sqr(g2d) + 0.75*g2u*traceYuAdjYu*Sqr(g2d
      ) - g2u*Lambdax*Sqr(g2d) + 0.15*g2u*Sqr(g1)*Sqr(g2d) + 4.25*g2u*Sqr(g2)*
      Sqr(g2d) - 4*g2d*gYd*gYu*Sqr(g2u) + 20*g2u*traceYuAdjYu*Sqr(g3) - 0.9375*
      Power(g2u,3)*Sqr(gYd) + 0.1875*g2u*Sqr(g1)*Sqr(gYd) + 0.9375*g2u*Sqr(g2)*
      Sqr(gYd) - 0.8125*g2u*Sqr(g2d)*Sqr(gYd) - 3.6875*Power(g2u,3)*Sqr(gYu) -
      1.125*g2u*traceYuAdjYu*Sqr(gYu) - g2u*Lambdax*Sqr(gYu) + 0.39375*g2u*Sqr(
      g1)*Sqr(gYu) + 3.46875*g2u*Sqr(g2)*Sqr(gYu) - 1.9375*g2u*Sqr(g2d)*Sqr(gYu
      ) - 1.5*g2u*Sqr(gYd)*Sqr(gYu) + 0.125*traceYdAdjYd*(-45*Power(g2u,3) - 24
      *g2d*gYd*gYu + 5*g2u*Sqr(g1) + 45*g2u*Sqr(g2) + 6*g2u*Sqr(g2d) + 160*g2u*
      Sqr(g3) - 9*g2u*Sqr(gYu)) + 0.125*traceYeAdjYe*(-15*Power(g2u,3) - 8*g2d*
      gYd*gYu + 15*g2u*Sqr(g1) + 15*g2u*Sqr(g2) + 2*g2u*Sqr(g2d) - 3*g2u*Sqr(
      gYu)) + 1.5*g2u*Sqr(Lambdax)));


   return beta_g2u;
}

/**
 * Calculates the three-loop beta function of g2u.
 *
 * @return three-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2u;

   beta_g2u = 0;


   return beta_g2u;
}

} // namespace flexiblesusy
