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

// File generated at Sun 26 Aug 2018 14:09:46

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of gYu.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_gYu;

   beta_gYu = Re(oneOver16PiSqr*(3*g2d*g2u*gYd + 1.5*gYu*Sqr(g2d) + 0.05*gYu*(-
      9*Sqr(g1) + 5*(12*traceYdAdjYd + 4*traceYeAdjYe + 12*traceYuAdjYu - 9*Sqr
      (g2) + 9*Sqr(g2u) + 8*Sqr(gYd) + 5*Sqr(gYu)))));


   return beta_gYu;
}

/**
 * Calculates the 2-loop beta function of gYu.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_gYu;

   beta_gYu = Re(0.00125*twoLoop*(468*gYu*Quad(g1) + 5*Sqr(g1)*(72*g2d*g2u*gYd
      + 90*gYu*Sqr(g2d) + gYu*(100*traceYdAdjYd + 300*traceYeAdjYe + 340*
      traceYuAdjYu - 216*Sqr(g2) + 189*Sqr(g2u) + 12*Sqr(gYd) + 309*Sqr(gYu)))
      - 25*(136*gYu*Quad(g2) - 3*Sqr(g2)*(136*g2d*g2u*gYd + 110*gYu*Sqr(g2d) +
      gYu*(60*traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu + 183*Sqr(g2u) +
      52*Sqr(gYd) + 55*Sqr(gYu))) + 2*(132*g2u*gYd*Cube(g2d) + 45*gYu*Quad(g2d)
      + 3*gYu*Sqr(g2d)*(14*Sqr(g2u) + 25*Sqr(gYd) + 9*Sqr(gYu)) + 12*g2d*g2u*
      gYd*(6*Sqr(g2u) + 3*Sqr(gYd) + 4*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu + Lambdax + 2*Sqr(gYu))) + gYu*(99*Quad(g2u) + 3*Sqr(g2u)*(2
      *(9*traceYdAdjYd + 3*traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax) + 25*Sqr(
      gYd) + 3*Sqr(gYu)) + 2*(18*Quad(gYd) + 6*Quad(gYu) + 3*(9*traceYdAdjYd +
      3*traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax)*Sqr(gYu) + 2*Sqr(gYd)*(21*
      traceYdAdjYd + 7*traceYeAdjYe + 21*traceYuAdjYu + 12*Lambdax + 15*Sqr(gYu
      )) - 2*(80*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3) - 3*(9*
      traceYdAdjYdYdAdjYd - 2*traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*
      traceYuAdjYuYuAdjYu - 2*Sqr(Lambdax)))))))));


   return beta_gYu;
}

/**
 * Calculates the 3-loop beta function of gYu.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYu;

   beta_gYu = 0;


   return beta_gYu;
}

/**
 * Calculates the 4-loop beta function of gYu.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYu;

   beta_gYu = 0;


   return beta_gYu;
}

} // namespace flexiblesusy
