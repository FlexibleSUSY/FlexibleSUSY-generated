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

// File generated at Tue 10 Oct 2017 21:09:50

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2d.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2d;

   beta_g2d = Re(oneOver16PiSqr*(g2u*gYd*gYu + 3*g2d*traceYdAdjYd + g2d*
      traceYeAdjYe + 3*g2d*traceYuAdjYu + 2.75*Cube(g2d) - 0.45*g2d*Sqr(g1) -
      8.25*g2d*Sqr(g2) + g2d*Sqr(g2u) + 0.75*g2d*Sqr(gYd) + 0.5*g2d*Sqr(gYu)));


   return beta_g2d;
}

/**
 * Calculates the 2-loop beta function of g2d.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2d;

   beta_g2d = Re(0.0004166666666666667*twoLoop*(1404*g2d*Quad(g1) + 15*
      Sqr(g1)*(24*g2u*gYd*gYu + 435*Cube(g2d) + 72*g2d*Sqr(g2) + g2d*(100*
      traceYdAdjYd + 300*traceYeAdjYe + 340*traceYuAdjYu + 24*Sqr(g2u) + 63*Sqr
      (gYd) + 30*Sqr(gYu))) - 25*(3272*g2d*Quad(g2) - 3*Sqr(g2)*(72*g2u*gYd*gYu
      + 875*Cube(g2d) + g2d*(136*Sqr(g2u) + 3*(37*Sqr(gYd) + 10*(6*
      traceYdAdjYd + 2*traceYeAdjYe + 6*traceYuAdjYu + Sqr(gYu))))) + 6*(56*
      Power5(g2d) + 64*g2u*gYd*gYu*Sqr(g2d) + 4*g2u*gYd*gYu*(12*traceYdAdjYd +
      4*traceYeAdjYe + 12*traceYuAdjYu + 4*Lambdax + 9*Sqr(g2u) + 6*Sqr(gYd) +
      5*Sqr(gYu)) + Cube(g2d)*(90*traceYdAdjYd + 30*traceYeAdjYe + 90*
      traceYuAdjYu + 80*Lambdax + 54*Sqr(g2u) + 59*Sqr(gYd) + 15*Sqr(gYu)) +
      g2d*(108*traceYdAdjYdYdAdjYd - 24*traceYdAdjYuYuAdjYd + 36*
      traceYeAdjYeYeAdjYe + 108*traceYuAdjYuYuAdjYu + 22*Quad(g2u) + 5*Quad(gYd
      ) + 9*Quad(gYu) - 320*traceYdAdjYd*Sqr(g3) - 320*traceYuAdjYu*Sqr(g3) + 2
      *Sqr(gYd)*(9*traceYdAdjYd + 3*traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax +
      12*Sqr(gYu)) + Sqr(g2u)*(-4*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 4*Lambdax) + 31*Sqr(gYd) + 13*Sqr(gYu)) - 24*Sqr(Lambdax))
      ))));


   return beta_g2d;
}

/**
 * Calculates the 3-loop beta function of g2d.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_g2d_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2d;

   beta_g2d = 0;


   return beta_g2d;
}

} // namespace flexiblesusy
