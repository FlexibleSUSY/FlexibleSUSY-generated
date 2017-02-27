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

// File generated at Mon 27 Feb 2017 13:24:36

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

   beta_gYd = Re(0.00125*twoLoop*(468*Power(g1,4)*gYd + 5*Sqr(g1)*(72*g2d
      *g2u*gYu - 216*gYd*Sqr(g2) + 189*gYd*Sqr(g2d) + gYd*(90*Sqr(g2u) + 309*
      Sqr(gYd) + 4*(25*traceYdAdjYd + 75*traceYeAdjYe + 85*traceYuAdjYu + 3*Sqr
      (gYu)))) - 25*(136*Power(g2,4)*gYd - 3*Sqr(g2)*(136*g2d*g2u*gYu + 183*gYd
      *Sqr(g2d) + gYd*(60*traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu +
      110*Sqr(g2u) + 55*Sqr(gYd) + 52*Sqr(gYu))) + 2*(99*Power(g2d,4)*gYd + 72*
      Power(g2d,3)*g2u*gYu + 12*g2d*g2u*gYu*(12*traceYdAdjYd + 4*traceYeAdjYe +
      12*traceYuAdjYu + 4*Lambdax + 11*Sqr(g2u) + 8*Sqr(gYd) + 3*Sqr(gYu)) + 3
      *gYd*Sqr(g2d)*(18*traceYdAdjYd + 6*traceYeAdjYe + 18*traceYuAdjYu + 16*
      Lambdax + 14*Sqr(g2u) + 3*Sqr(gYd) + 25*Sqr(gYu)) + gYd*(45*Power(g2u,4)
      + 3*Sqr(g2u)*(9*Sqr(gYd) + 25*Sqr(gYu)) + 2*(6*Power(gYd,4) + 3*Sqr(gYd)*
      (9*traceYdAdjYd + 3*traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax + 10*Sqr(
      gYu)) + 2*(9*Power(gYu,4) - 80*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + (
      21*traceYdAdjYd + 7*traceYeAdjYe + 21*traceYuAdjYu + 12*Lambdax)*Sqr(gYu)
      + 3*(9*traceYdAdjYdYdAdjYd - 2*traceYdAdjYuYuAdjYd + 3*
      traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2*Sqr(Lambdax)))))))));


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
