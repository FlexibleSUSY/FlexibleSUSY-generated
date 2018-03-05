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

// File generated at Mon 5 Mar 2018 16:34:43

#include "E6SSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vd.
 *
 * @return 1-loop beta function
 */
double E6SSMEFTHiggs_susy_parameters::calc_beta_vd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_vd;

   beta_vd = Re(0.05*oneOver16PiSqr*vd*(-60*traceYdAdjYd - 20*
      traceYeAdjYe - 20*AbsSqr(Lambdax) + 6*Sqr(g1) + 30*Sqr(g2) + 9*Sqr(gN)));


   return beta_vd;
}

/**
 * Calculates the 2-loop beta function of vd.
 *
 * @return 2-loop beta function
 */
double E6SSMEFTHiggs_susy_parameters::calc_beta_vd_2_loop(const Susy_traces& susy_traces) const
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


   double beta_vd;

   beta_vd = Re(-0.00125*twoLoop*vd*(-7200*traceYdAdjYdYdAdjYd - 2400*
      traceYdAdjYuYuAdjYd - 2400*traceYeAdjYeYeAdjYe + 1188*Quad(g1) + 2900*
      Quad(g2) + 1773*Quad(gN) + 3600*traceYdAdjYd*Sqr(g2) + 1200*traceYeAdjYe*
      Sqr(g2) + 12800*traceYdAdjYd*Sqr(g3) + 4*Sqr(g1)*(100*(traceYdAdjYd + 3*
      traceYeAdjYe) + 90*Sqr(g2) - 9*Sqr(gN)) + 600*traceYdAdjYd*Sqr(gN) + 200*
      traceYeAdjYe*Sqr(gN) + 540*Sqr(g2)*Sqr(gN) + 40*AbsSqr(Lambdax)*(-60*
      traceKappaAdjKappa - 40*traceLambda12AdjLambda12 - 60*traceYuAdjYu + 6*
      Sqr(g1) + 30*Sqr(g2) + 29*Sqr(gN)) - 2400*Sqr(Conj(Lambdax))*Sqr(Lambdax)
      ));


   return beta_vd;
}

/**
 * Calculates the 3-loop beta function of vd.
 *
 * @return 3-loop beta function
 */
double E6SSMEFTHiggs_susy_parameters::calc_beta_vd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vd;

   beta_vd = 0;


   return beta_vd;
}

/**
 * Calculates the 4-loop beta function of vd.
 *
 * @return 4-loop beta function
 */
double E6SSMEFTHiggs_susy_parameters::calc_beta_vd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vd;

   beta_vd = 0;


   return beta_vd;
}

} // namespace flexiblesusy
