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

// File generated at Tue 8 Sep 2015 12:36:11

#include "E6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of gN.
 *
 * @return one-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_gN_one_loop(const Susy_traces& susy_traces) const
{


   double beta_gN;

   beta_gN = Re(9.4*Power(gN,3)*oneOver16PiSqr);


   return beta_gN;
}

/**
 * Calculates the two-loop beta function of gN.
 *
 * @return two-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_gN_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   double beta_gN;

   beta_gN = Re(0.02*Power(gN,3)*twoLoop*(-285*traceKappaAdjKappa - 190*
      traceLambda12AdjLambda12 - 210*traceYdAdjYd - 70*traceYeAdjYe - 90*
      traceYuAdjYu - 190*AbsSqr(Lambdax) + 162*Sqr(g1) + 510*Sqr(g2) + 1200*Sqr
      (g3) + 458*Sqr(gN)));


   return beta_gN;
}

/**
 * Calculates the three-loop beta function of gN.
 *
 * @return three-loop beta function
 */
double E6SSM_susy_parameters::calc_beta_gN_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gN;

   beta_gN = 0;


   return beta_gN;
}

} // namespace flexiblesusy
