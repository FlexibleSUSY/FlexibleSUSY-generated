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

// File generated at Mon 5 Mar 2018 17:41:58

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambdax.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_Lambdax_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(oneOver16PiSqr*(-4*g2d*g2u*gYd*gYu - 12*
      traceYdAdjYdYdAdjYd - 4*traceYeAdjYeYeAdjYe - 12*traceYuAdjYuYuAdjYu + 12
      *traceYdAdjYd*Lambdax + 4*traceYeAdjYe*Lambdax + 12*traceYuAdjYu*Lambdax
      + 0.27*Quad(g1) + 2.25*Quad(g2) - 5*Quad(g2d) - 5*Quad(g2u) - Quad(gYd) -
      Quad(gYu) - 9*Lambdax*Sqr(g2) + 0.9*Sqr(g1)*(-2*Lambdax + Sqr(g2)) + 6*
      Lambdax*Sqr(g2d) + 6*Lambdax*Sqr(g2u) - 2*Sqr(g2d)*Sqr(g2u) + 2*Lambdax*
      Sqr(gYd) - 2*Sqr(g2d)*Sqr(gYd) + 2*Lambdax*Sqr(gYu) - 2*Sqr(g2u)*Sqr(gYu)
      - 2*Sqr(gYd)*Sqr(gYu) + 12*Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the 2-loop beta function of Lambdax.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_Lambdax_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYdYdAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe =
      TRACE_STRUCT.traceYeAdjYeYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYuYuAdjYu =
      TRACE_STRUCT.traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_Lambdax;

   const double beta_Lambdax_1 = Re(-0.001*twoLoop*(3699*Power6(g1) + 45*
      Quad(g1)*(-20*traceYdAdjYd + 100*traceYeAdjYe + 76*traceYuAdjYu - 223*
      Lambdax + 197*Sqr(g2) + 6*Sqr(g2d) + 6*Sqr(g2u) + 2*Sqr(gYd) + 2*Sqr(gYu)
      ) + 25*Sqr(g1)*(385*Quad(g2) - 6*Sqr(g2)*(36*traceYdAdjYd + 44*
      traceYeAdjYe + 84*traceYuAdjYu + 39*Lambdax + 42*Sqr(g2d) + 42*Sqr(g2u) -
      2*Sqr(gYd) - 2*Sqr(gYu)) - 2*(32*traceYdAdjYdYdAdjYd - 96*
      traceYeAdjYeYeAdjYe - 64*traceYuAdjYuYuAdjYu + 50*traceYdAdjYd*Lambdax +
      150*traceYeAdjYe*Lambdax + 170*traceYuAdjYu*Lambdax + 45*Lambdax*Sqr(g2d)
      + 45*Lambdax*Sqr(g2u) + 15*Lambdax*Sqr(gYd) + 15*Lambdax*Sqr(gYu) + 216*
      Sqr(Lambdax))) - 125*(209*Power6(g2) - Quad(g2)*(36*traceYdAdjYd + 12*
      traceYeAdjYe + 36*traceYuAdjYu - 47*Lambdax + 306*Sqr(g2d) + 306*Sqr(g2u)
      + 6*Sqr(gYd) + 6*Sqr(gYu)) - 2*Sqr(g2)*(32*g2d*g2u*gYd*gYu + 80*Quad(g2d
      ) + 80*Quad(g2u) + Sqr(g2d)*(-165*Lambdax + 32*Sqr(g2u) + 16*Sqr(gYd)) -
      3*Lambdax*(30*traceYdAdjYd + 10*traceYeAdjYe + 30*traceYuAdjYu + 72*
      Lambdax + 5*Sqr(gYd) + 5*Sqr(gYu)) + Sqr(g2u)*(-165*Lambdax + 16*Sqr(gYu)
      )) + 2*(240*traceYdAdjYdYdAdjYdYdAdjYd - 96*traceYdAdjYdYdAdjYuYuAdjYd +
      48*traceYdAdjYuYuAdjYdYdAdjYd - 48*traceYdAdjYuYuAdjYuYuAdjYd + 80*
      traceYeAdjYeYeAdjYeYeAdjYe + 240*traceYuAdjYuYuAdjYuYuAdjYu + 76*g2u*gYd*
      gYu*Cube(g2d) - 12*traceYdAdjYdYdAdjYd*Lambdax - 168*traceYdAdjYuYuAdjYd*
      Lambdax - 4*traceYeAdjYeYeAdjYe*Lambdax - 12*traceYuAdjYuYuAdjYu*Lambdax
      + 94*Power6(g2d) + 94*Power6(g2u) + 10*Power6(gYd) + 10*Power6(gYu) -
      Lambdax*Quad(gYd) - Lambdax*Quad(gYu) - 256*traceYdAdjYdYdAdjYd*Sqr(g3) -
      256*traceYuAdjYuYuAdjYu*Sqr(g3) + 320*traceYdAdjYd*Lambdax*Sqr(g3) + 320
      *traceYuAdjYu*Lambdax*Sqr(g3) + 34*Quad(gYu)*Sqr(gYd) + Quad(g2d)*(-5*
      Lambdax + 14*Sqr(g2u) + 22*Sqr(gYd)) + 34*Quad(gYd)*Sqr(gYu) + 12*Lambdax
      *Sqr(gYd)*Sqr(gYu) + 4*g2d*g2u*gYd*gYu*(20*Lambdax + 19*Sqr(g2u) + 21*Sqr
      (gYd) + 21*Sqr(gYu)) + Quad(g2u)*(-5*Lambdax + 22*Sqr(gYu)) + 2*Sqr(g2u)*
      (17*Quad(gYu) - Lambdax*Sqr(gYu) + 19*Sqr(gYd)*Sqr(gYu) - 72*Sqr(Lambdax)
      ) + 2*Sqr(g2d)*(7*Quad(g2u) + 17*Quad(gYd) + Sqr(gYd)*(-Lambdax + 19*Sqr(
      gYu)) + Sqr(g2u)*(-22*Lambdax + 21*Sqr(gYd) + 21*Sqr(gYu)) - 72*Sqr(
      Lambdax)) - 48*Sqr(gYd)*Sqr(Lambdax) - 48*Sqr(gYu)*Sqr(Lambdax)))));
   const double beta_Lambdax_2 = Re(-6*twoLoop*(12*traceYdAdjYd + 4*
      traceYeAdjYe + 12*traceYuAdjYu + 13*Lambdax)*Sqr(Lambdax));

   beta_Lambdax = beta_Lambdax_1 + beta_Lambdax_2;


   return beta_Lambdax;
}

/**
 * Calculates the 3-loop beta function of Lambdax.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_Lambdax_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

/**
 * Calculates the 4-loop beta function of Lambdax.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_Lambdax_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
