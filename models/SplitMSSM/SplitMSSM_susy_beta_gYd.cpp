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

// File generated at Sun 26 Aug 2018 14:09:43

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of gYd.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_gYd;

   beta_gYd = Re(oneOver16PiSqr*(3*g2d*g2u*gYu + 3*gYd*traceYdAdjYd + gYd*
      traceYeAdjYe + 3*gYd*traceYuAdjYu + 1.25*Cube(gYd) - 0.45*gYd*Sqr(g1) -
      2.25*gYd*Sqr(g2) + 2.25*gYd*Sqr(g2d) + 1.5*gYd*Sqr(g2u) + 2*gYd*Sqr(gYu))
      );


   return beta_gYd;
}

/**
 * Calculates the 2-loop beta function of gYd.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_gYd;

   beta_gYd = Re(0.00125*twoLoop*(468*gYd*Quad(g1) + 5*Sqr(g1)*(72*g2d*g2u*gYu
      - 216*gYd*Sqr(g2) + 189*gYd*Sqr(g2d) + gYd*(90*Sqr(g2u) + 309*Sqr(gYd) +
      4*(25*traceYdAdjYd + 75*traceYeAdjYe + 85*traceYuAdjYu + 3*Sqr(gYu)))) -
      25*(136*gYd*Quad(g2) - 3*Sqr(g2)*(136*g2d*g2u*gYu + 183*gYd*Sqr(g2d) +
      gYd*(60*traceYdAdjYd + 20*traceYeAdjYe + 60*traceYuAdjYu + 110*Sqr(g2u) +
      55*Sqr(gYd) + 52*Sqr(gYu))) + 2*(72*g2u*gYu*Cube(g2d) + 99*gYd*Quad(g2d)
      + 12*g2d*g2u*gYu*(12*traceYdAdjYd + 4*traceYeAdjYe + 12*traceYuAdjYu + 4*
      Lambdax + 11*Sqr(g2u) + 8*Sqr(gYd) + 3*Sqr(gYu)) + 3*gYd*Sqr(g2d)*(18*
      traceYdAdjYd + 6*traceYeAdjYe + 18*traceYuAdjYu + 16*Lambdax + 14*Sqr(g2u
      ) + 3*Sqr(gYd) + 25*Sqr(gYu)) + gYd*(45*Quad(g2u) + 3*Sqr(g2u)*(9*Sqr(gYd
      ) + 25*Sqr(gYu)) + 2*(6*Quad(gYd) + 3*Sqr(gYd)*(9*traceYdAdjYd + 3*
      traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax + 10*Sqr(gYu)) + 2*(9*Quad(gYu)
      - 80*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + (21*traceYdAdjYd + 7*
      traceYeAdjYe + 21*traceYuAdjYu + 12*Lambdax)*Sqr(gYu) + 3*(9*
      traceYdAdjYdYdAdjYd - 2*traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*
      traceYuAdjYuYuAdjYu - 2*Sqr(Lambdax)))))))));


   return beta_gYd;
}

/**
 * Calculates the 3-loop beta function of gYd.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYd;

   beta_gYd = 0;


   return beta_gYd;
}

/**
 * Calculates the 4-loop beta function of gYd.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYd;

   beta_gYd = 0;


   return beta_gYd;
}

} // namespace flexiblesusy
