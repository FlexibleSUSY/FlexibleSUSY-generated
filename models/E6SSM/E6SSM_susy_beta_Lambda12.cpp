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

// File generated at Fri 20 Oct 2017 08:51:15

#include "E6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda12.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSM_susy_parameters::calc_beta_Lambda12_1_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,2,2> beta_Lambda12;

   beta_Lambda12 = (oneOver16PiSqr*(Lambda12*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 2*AbsSqr(Lambdax) - 0.6*Sqr(g1) - 3*Sqr(g2) -
      1.9*Sqr(gN)) + 2*(Lambda12*(Lambda12).adjoint()*Lambda12))).real();


   return beta_Lambda12;
}

/**
 * Calculates the 2-loop beta function of Lambda12.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSM_susy_parameters::calc_beta_Lambda12_2_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,2,2> beta_Lambda12;

   beta_Lambda12 = (twoLoop*(0.005*Lambda12*(-1200*
      traceKappaAdjKappaKappaAdjKappa - 800*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 1188*Quad(g1) + 3300*Quad(
      g2) + 3933*Quad(gN) + 1200*traceLambda12AdjLambda12*Sqr(g2) + 3200*
      traceKappaAdjKappa*Sqr(g3) + 80*AbsSqr(Lambdax)*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 15*traceYuAdjYu + 3*Sqr(g1) + 15*Sqr(g2) - 3*Sqr(gN)) -
      360*traceKappaAdjKappa*Sqr(gN) - 240*traceLambda12AdjLambda12*Sqr(gN) +
      390*Sqr(g2)*Sqr(gN) + 2*Sqr(g1)*(80*traceKappaAdjKappa + 120*
      traceLambda12AdjLambda12 + 180*Sqr(g2) + 27*Sqr(gN)) - 800*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)) + (-6*traceKappaAdjKappa - 4*
      traceLambda12AdjLambda12 - 4*AbsSqr(Lambdax) + 2.5*Sqr(gN))*(Lambda12*(
      Lambda12).adjoint()*Lambda12) - 2*(Lambda12*(Lambda12).adjoint()*Lambda12
      *(Lambda12).adjoint()*Lambda12))).real();


   return beta_Lambda12;
}

/**
 * Calculates the 3-loop beta function of Lambda12.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSM_susy_parameters::calc_beta_Lambda12_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,2,2> beta_Lambda12;

   beta_Lambda12 = ZEROMATRIX(2,2);


   return beta_Lambda12;
}

} // namespace flexiblesusy
