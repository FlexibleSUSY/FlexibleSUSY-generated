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

// File generated at Fri 20 Oct 2017 08:51:16

#include "E6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambdax.
 *
 * @return 1-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_Lambdax_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_Lambdax;

   beta_Lambdax = Re(0.1*oneOver16PiSqr*Lambdax*(30*traceKappaAdjKappa +
      20*traceLambda12AdjLambda12 + 30*traceYdAdjYd + 10*traceYeAdjYe + 30*
      traceYuAdjYu + 40*AbsSqr(Lambdax) - 6*Sqr(g1) - 30*Sqr(g2) - 19*Sqr(gN)))
      ;


   return beta_Lambdax;
}

/**
 * Calculates the 2-loop beta function of Lambdax.
 *
 * @return 2-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_Lambdax_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   double beta_Lambdax;

   beta_Lambdax = Re(-0.005*twoLoop*Lambdax*(1200*
      traceKappaAdjKappaKappaAdjKappa + 800*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 1800*traceYdAdjYdYdAdjYd +
      1200*traceYdAdjYuYuAdjYd + 600*traceYeAdjYeYeAdjYe + 1800*
      traceYuAdjYuYuAdjYu - 1188*Quad(g1) - 3300*Quad(g2) - 3933*Quad(gN) -
      1200*traceLambda12AdjLambda12*Sqr(g2) - 3200*traceKappaAdjKappa*Sqr(g3) -
      3200*traceYdAdjYd*Sqr(g3) - 3200*traceYuAdjYu*Sqr(g3) + 360*
      traceKappaAdjKappa*Sqr(gN) + 240*traceLambda12AdjLambda12*Sqr(gN) + 120*
      traceYdAdjYd*Sqr(gN) + 40*traceYeAdjYe*Sqr(gN) + 60*traceYuAdjYu*Sqr(gN)
      - 390*Sqr(g2)*Sqr(gN) - 20*AbsSqr(Lambdax)*(-60*traceKappaAdjKappa - 40*
      traceLambda12AdjLambda12 - 90*traceYdAdjYd - 30*traceYeAdjYe - 90*
      traceYuAdjYu + 12*Sqr(g1) + 60*Sqr(g2) + 13*Sqr(gN)) - 2*Sqr(g1)*(40*(2*
      traceKappaAdjKappa + 3*traceLambda12AdjLambda12 - traceYdAdjYd + 3*
      traceYeAdjYe + 2*traceYuAdjYu) + 180*Sqr(g2) + 27*Sqr(gN)) + 2000*Sqr(
      Conj(Lambdax))*Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the 3-loop beta function of Lambdax.
 *
 * @return 3-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_Lambdax_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
