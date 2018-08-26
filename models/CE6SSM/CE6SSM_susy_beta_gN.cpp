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

// File generated at Sun 26 Aug 2018 13:46:25

#include "CE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of gN.
 *
 * @return 1-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_gN_1_loop(const Susy_traces& susy_traces) const
{


   double beta_gN;

   beta_gN = Re(9.4*oneOver16PiSqr*Cube(gN));


   return beta_gN;
}

/**
 * Calculates the 2-loop beta function of gN.
 *
 * @return 2-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_gN_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;


   double beta_gN;

   beta_gN = Re(0.02*twoLoop*Cube(gN)*(-285*traceKappaAdjKappa - 190*
      traceLambda12AdjLambda12 - 210*traceYdAdjYd - 70*traceYeAdjYe - 90*
      traceYuAdjYu - 190*AbsSqr(Lambdax) + 162*Sqr(g1) + 510*Sqr(g2) + 1200*Sqr
      (g3) + 458*Sqr(gN)));


   return beta_gN;
}

/**
 * Calculates the 3-loop beta function of gN.
 *
 * @return 3-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_gN_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gN;

   beta_gN = 0;


   return beta_gN;
}

/**
 * Calculates the 4-loop beta function of gN.
 *
 * @return 4-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_gN_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gN;

   beta_gN = 0;


   return beta_gN;
}

} // namespace flexiblesusy
