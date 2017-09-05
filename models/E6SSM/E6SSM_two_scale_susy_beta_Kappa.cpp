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

// File generated at Tue 5 Sep 2017 11:07:41

#include "E6SSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Kappa.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Kappa_one_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = (oneOver16PiSqr*(Kappa*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 2*AbsSqr(Lambdax) - 0.26666666666666666*Sqr(g1
      ) - 5.333333333333333*Sqr(g3) - 1.9*Sqr(gN)) + 2*(Kappa*(Kappa).adjoint()
      *Kappa))).real();


   return beta_Kappa;
}

/**
 * Calculates the two-loop beta function of Kappa.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Kappa_two_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = (twoLoop*(Kappa*(2.5955555555555554*Power(g1,4) +
      14.222222222222221*Power(g3,4) + 19.665*Power(gN,4) - 6*
      traceKappaAdjKappaKappaAdjKappa - 4*
      traceLambda12AdjLambda12Lambda12AdjLambda12 + 6*traceLambda12AdjLambda12*
      Sqr(g2) + 16*traceKappaAdjKappa*Sqr(g3) + 0.4*AbsSqr(Lambdax)*(-15*
      traceYdAdjYd - 5*traceYeAdjYe - 15*traceYuAdjYu + 3*Sqr(g1) + 15*Sqr(g2)
      - 3*Sqr(gN)) - 1.8*traceKappaAdjKappa*Sqr(gN) - 1.2*
      traceLambda12AdjLambda12*Sqr(gN) + 3.466666666666667*Sqr(g3)*Sqr(gN) +
      0.0044444444444444444*Sqr(g1)*(180*traceKappaAdjKappa + 270*
      traceLambda12AdjLambda12 + 320*Sqr(g3) + 57*Sqr(gN)) - 4*Sqr(Conj(Lambdax
      ))*Sqr(Lambdax)) + (-6*traceKappaAdjKappa - 4*traceLambda12AdjLambda12 -
      4*AbsSqr(Lambdax) + 2.5*Sqr(gN))*(Kappa*(Kappa).adjoint()*Kappa) - 2*(
      Kappa*(Kappa).adjoint()*Kappa*(Kappa).adjoint()*Kappa))).real();


   return beta_Kappa;
}

/**
 * Calculates the three-loop beta function of Kappa.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Kappa_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = ZEROMATRIX(3,3);


   return beta_Kappa;
}

} // namespace flexiblesusy
