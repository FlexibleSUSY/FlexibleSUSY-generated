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

// File generated at Sun 4 Aug 2019 17:10:41

#include "CE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambdax.
 *
 * @return 1-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_Lambdax_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;


   double beta_Lambdax;

   beta_Lambdax = Re(0.1*oneOver16PiSqr*Lambdax*(30*traceKappaAdjKappa + 20*
      traceLambda12AdjLambda12 + 30*traceYdAdjYd + 10*traceYeAdjYe + 30*
      traceYuAdjYu + 40*AbsSqr(Lambdax) - 6*Sqr(g1) - 30*Sqr(g2) - 19*Sqr(gN)))
      ;


   return beta_Lambdax;
}

/**
 * Calculates the 2-loop beta function of Lambdax.
 *
 * @return 2-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_Lambdax_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceKappaAdjKappaKappaAdjKappa = TRACE_STRUCT.
      traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12Lambda12AdjLambda12;


   double beta_Lambdax;

   beta_Lambdax = Re(-0.005*twoLoop*Lambdax*(1200*
      traceKappaAdjKappaKappaAdjKappa + 800*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 1800*traceYdAdjYdYdAdjYd +
      1200*traceYdAdjYuYuAdjYd + 600*traceYeAdjYeYeAdjYe + 1800*
      traceYuAdjYuYuAdjYu + 1200*traceKappaAdjKappa*AbsSqr(Lambdax) + 800*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) + 1800*traceYdAdjYd*AbsSqr(
      Lambdax) + 600*traceYeAdjYe*AbsSqr(Lambdax) + 1800*traceYuAdjYu*AbsSqr(
      Lambdax) - 1188*Quad(g1) - 3300*Quad(g2) - 3933*Quad(gN) - 160*
      traceKappaAdjKappa*Sqr(g1) - 240*traceLambda12AdjLambda12*Sqr(g1) + 80*
      traceYdAdjYd*Sqr(g1) - 240*traceYeAdjYe*Sqr(g1) - 160*traceYuAdjYu*Sqr(g1
      ) - 240*AbsSqr(Lambdax)*Sqr(g1) - 1200*traceLambda12AdjLambda12*Sqr(g2) -
      1200*AbsSqr(Lambdax)*Sqr(g2) - 360*Sqr(g1)*Sqr(g2) - 3200*
      traceKappaAdjKappa*Sqr(g3) - 3200*traceYdAdjYd*Sqr(g3) - 3200*
      traceYuAdjYu*Sqr(g3) + 360*traceKappaAdjKappa*Sqr(gN) + 240*
      traceLambda12AdjLambda12*Sqr(gN) + 120*traceYdAdjYd*Sqr(gN) + 40*
      traceYeAdjYe*Sqr(gN) + 60*traceYuAdjYu*Sqr(gN) - 260*AbsSqr(Lambdax)*Sqr(
      gN) - 54*Sqr(g1)*Sqr(gN) - 390*Sqr(g2)*Sqr(gN) + 2000*Sqr(Conj(Lambdax))*
      Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the 3-loop beta function of Lambdax.
 *
 * @return 3-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_Lambdax_3_loop(const Susy_traces& susy_traces) const
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
double CE6SSM_susy_parameters::calc_beta_Lambdax_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

/**
 * Calculates the 5-loop beta function of Lambdax.
 *
 * @return 5-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_Lambdax_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
