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

// File generated at Mon 5 Mar 2018 17:42:04

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2u.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(oneOver16PiSqr*(g2d*gYd*gYu + 3*g2u*traceYdAdjYd + g2u*
      traceYeAdjYe + 3*g2u*traceYuAdjYu + 2.75*Cube(g2u) - 0.45*g2u*Sqr(g1) -
      8.25*g2u*Sqr(g2) + g2u*Sqr(g2d) + 0.5*g2u*Sqr(gYd) + 0.75*g2u*Sqr(gYu)));


   return beta_g2u;
}

/**
 * Calculates the 2-loop beta function of g2u.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2u;

   beta_g2u = Re(0.0004166666666666667*twoLoop*(1404*g2u*Quad(g1) + 15*
      Sqr(g1)*(24*g2d*gYd*gYu + 72*g2u*Sqr(g2) + 24*g2u*Sqr(g2d) + g2u*(100*
      traceYdAdjYd + 300*traceYeAdjYe + 340*traceYuAdjYu + 435*Sqr(g2u) + 30*
      Sqr(gYd) + 63*Sqr(gYu))) - 25*(3272*g2u*Quad(g2) - 3*Sqr(g2)*(72*g2d*gYd*
      gYu + 136*g2u*Sqr(g2d) + g2u*(180*traceYdAdjYd + 60*traceYeAdjYe + 180*
      traceYuAdjYu + 875*Sqr(g2u) + 30*Sqr(gYd) + 111*Sqr(gYu))) + 6*(36*gYd*
      gYu*Cube(g2d) + 22*g2u*Quad(g2d) + 4*g2d*gYd*gYu*(12*traceYdAdjYd + 4*
      traceYeAdjYe + 12*traceYuAdjYu + 4*Lambdax + 16*Sqr(g2u) + 5*Sqr(gYd) + 6
      *Sqr(gYu)) + g2u*Sqr(g2d)*(-12*traceYdAdjYd - 4*traceYeAdjYe - 12*
      traceYuAdjYu + 16*Lambdax + 54*Sqr(g2u) + 13*Sqr(gYd) + 31*Sqr(gYu)) +
      g2u*(108*traceYdAdjYdYdAdjYd - 24*traceYdAdjYuYuAdjYd + 36*
      traceYeAdjYeYeAdjYe + 108*traceYuAdjYuYuAdjYu + 56*Quad(g2u) + 9*Quad(gYd
      ) + 5*Quad(gYu) - 320*traceYdAdjYd*Sqr(g3) - 320*traceYuAdjYu*Sqr(g3) +
      18*traceYdAdjYd*Sqr(gYu) + 6*traceYeAdjYe*Sqr(gYu) + 18*traceYuAdjYu*Sqr(
      gYu) + 16*Lambdax*Sqr(gYu) + 24*Sqr(gYd)*Sqr(gYu) + Sqr(g2u)*(90*
      traceYdAdjYd + 30*traceYeAdjYe + 90*traceYuAdjYu + 80*Lambdax + 15*Sqr(
      gYd) + 59*Sqr(gYu)) - 24*Sqr(Lambdax))))));


   return beta_g2u;
}

/**
 * Calculates the 3-loop beta function of g2u.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2u;

   beta_g2u = 0;


   return beta_g2u;
}

/**
 * Calculates the 4-loop beta function of g2u.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2u_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2u;

   beta_g2u = 0;


   return beta_g2u;
}

} // namespace flexiblesusy
