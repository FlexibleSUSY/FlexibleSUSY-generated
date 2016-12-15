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

// File generated at Thu 15 Dec 2016 12:42:28

#include "HSSUSY_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambdax.
 *
 * @return one-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(oneOver16PiSqr*(0.27*Power(g1,4) + 2.25*Power(g2,4)
      - 9*Lambdax*Sqr(g2) + 0.9*Sqr(g1)*(-2*Lambdax + Sqr(g2)) + 4*(-3*
      traceYdAdjYdYdAdjYd - traceYeAdjYeYeAdjYe - 3*traceYuAdjYuYuAdjYu + 3*
      traceYdAdjYd*Lambdax + traceYeAdjYe*Lambdax + 3*traceYuAdjYu*Lambdax + 3*
      Sqr(Lambdax))));


   return beta_Lambdax;
}

/**
 * Calculates the two-loop beta function of Lambdax.
 *
 * @return two-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_two_loop(const Susy_traces& susy_traces) const
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

   beta_Lambdax = Re(twoLoop*(-3.411*Power(g1,6) + 38.125*Power(g2,6) +
      60*traceYdAdjYdYdAdjYdYdAdjYd - 24*traceYdAdjYdYdAdjYuYuAdjYd + 12*
      traceYdAdjYuYuAdjYdYdAdjYd - 12*traceYdAdjYuYuAdjYuYuAdjYd + 20*
      traceYeAdjYeYeAdjYeYeAdjYe + 60*traceYuAdjYuYuAdjYuYuAdjYu - 3*
      traceYdAdjYdYdAdjYd*Lambdax - 42*traceYdAdjYuYuAdjYd*Lambdax -
      traceYeAdjYeYeAdjYe*Lambdax - 3*traceYuAdjYuYuAdjYu*Lambdax - 78*Power(
      Lambdax,3) - 0.125*Power(g2,4)*(36*traceYdAdjYd + 12*traceYeAdjYe + 36*
      traceYuAdjYu + 73*Lambdax) + 1.5*Lambdax*(15*traceYdAdjYd + 5*
      traceYeAdjYe + 15*traceYuAdjYu + 36*Lambdax)*Sqr(g2) - 0.015*Power(g1,4)*
      (-60*traceYdAdjYd + 300*traceYeAdjYe + 228*traceYuAdjYu - 629*Lambdax +
      559*Sqr(g2)) - 64*traceYdAdjYdYdAdjYd*Sqr(g3) - 64*traceYuAdjYuYuAdjYu*
      Sqr(g3) + 80*traceYdAdjYd*Lambdax*Sqr(g3) + 80*traceYuAdjYu*Lambdax*Sqr(
      g3) - 72*traceYdAdjYd*Sqr(Lambdax) - 24*traceYeAdjYe*Sqr(Lambdax) - 72*
      traceYuAdjYu*Sqr(Lambdax) - 0.025*Sqr(g1)*(289*Power(g2,4) - 6*(36*
      traceYdAdjYd + 44*traceYeAdjYe + 84*traceYuAdjYu + 39*Lambdax)*Sqr(g2) -
      4*(16*traceYdAdjYdYdAdjYd - 48*traceYeAdjYeYeAdjYe - 32*
      traceYuAdjYuYuAdjYu + 25*traceYdAdjYd*Lambdax + 75*traceYeAdjYe*Lambdax +
      85*traceYuAdjYu*Lambdax + 108*Sqr(Lambdax)))));


   return beta_Lambdax;
}

/**
 * Calculates the three-loop beta function of Lambdax.
 *
 * @return three-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
