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

// File generated at Mon 5 Mar 2018 16:34:44

#include "E6SSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vs.
 *
 * @return 1-loop beta function
 */
double E6SSMEFTHiggs_susy_parameters::calc_beta_vs_1_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_vs;

   beta_vs = Re(0.25*oneOver16PiSqr*vs*(-12*traceKappaAdjKappa - 8*
      traceLambda12AdjLambda12 - 8*AbsSqr(Lambdax) + 5*Sqr(gN)));


   return beta_vs;
}

/**
 * Calculates the 2-loop beta function of vs.
 *
 * @return 2-loop beta function
 */
double E6SSMEFTHiggs_susy_parameters::calc_beta_vs_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;


   double beta_vs;

   beta_vs = Re(-0.00625*twoLoop*vs*(1065*Quad(gN) + 64*((2*
      traceKappaAdjKappa + 3*traceLambda12AdjLambda12)*Sqr(g1) + 5*(-3*
      traceKappaAdjKappaKappaAdjKappa - 2*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 3*traceLambda12AdjLambda12*
      Sqr(g2) + 8*traceKappaAdjKappa*Sqr(g3))) + 104*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12)*Sqr(gN) + 16*AbsSqr(Lambdax)*(-60*traceYdAdjYd
      - 20*traceYeAdjYe - 60*traceYuAdjYu + 12*Sqr(g1) + 60*Sqr(g2) + 13*Sqr(gN
      )) - 640*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_vs;
}

/**
 * Calculates the 3-loop beta function of vs.
 *
 * @return 3-loop beta function
 */
double E6SSMEFTHiggs_susy_parameters::calc_beta_vs_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vs;

   beta_vs = 0;


   return beta_vs;
}

/**
 * Calculates the 4-loop beta function of vs.
 *
 * @return 4-loop beta function
 */
double E6SSMEFTHiggs_susy_parameters::calc_beta_vs_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vs;

   beta_vs = 0;


   return beta_vs;
}

} // namespace flexiblesusy
