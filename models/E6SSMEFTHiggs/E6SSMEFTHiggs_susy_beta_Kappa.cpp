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

// File generated at Fri 10 Apr 2020 18:30:21

#include "E6SSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Kappa.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Kappa_1_loop(const Susy_traces& susy_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = (oneOver16PiSqr*(-0.03333333333333333*Kappa*(-90*
      traceKappaAdjKappa - 60*traceLambda12AdjLambda12 - 60*AbsSqr(Lambdax) + 8
      *Sqr(g1) + 160*Sqr(g3) + 57*Sqr(gN)) + 2*(Kappa*(Kappa).adjoint()*Kappa))
      ).real();


   return beta_Kappa;
}

/**
 * Calculates the 2-loop beta function of Kappa.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Kappa_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;
   const double traceKappaAdjKappaKappaAdjKappa = TRACE_STRUCT.
      traceKappaAdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12Lambda12AdjLambda12;


   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = (twoLoop*(0.0005555555555555556*Kappa*(-10800*
      traceKappaAdjKappaKappaAdjKappa - 7200*
      traceLambda12AdjLambda12Lambda12AdjLambda12 - 10800*traceYdAdjYd*AbsSqr(
      Lambdax) - 3600*traceYeAdjYe*AbsSqr(Lambdax) - 10800*traceYuAdjYu*AbsSqr(
      Lambdax) + 4672*Quad(g1) + 25600*Quad(g3) + 35397*Quad(gN) + 1440*
      traceKappaAdjKappa*Sqr(g1) + 2160*traceLambda12AdjLambda12*Sqr(g1) + 2160
      *AbsSqr(Lambdax)*Sqr(g1) + 10800*traceLambda12AdjLambda12*Sqr(g2) + 10800
      *AbsSqr(Lambdax)*Sqr(g2) + 28800*traceKappaAdjKappa*Sqr(g3) + 2560*Sqr(g1
      )*Sqr(g3) - 3240*traceKappaAdjKappa*Sqr(gN) - 2160*
      traceLambda12AdjLambda12*Sqr(gN) - 2160*AbsSqr(Lambdax)*Sqr(gN) + 456*Sqr
      (g1)*Sqr(gN) + 6240*Sqr(g3)*Sqr(gN) - 7200*Sqr(Conj(Lambdax))*Sqr(Lambdax
      )) + 0.5*(-12*traceKappaAdjKappa - 8*traceLambda12AdjLambda12 - 8*AbsSqr(
      Lambdax) + 5*Sqr(gN))*(Kappa*(Kappa).adjoint()*Kappa) - 2*(Kappa*(Kappa).
      adjoint()*Kappa*(Kappa).adjoint()*Kappa))).real();


   return beta_Kappa;
}

/**
 * Calculates the 3-loop beta function of Kappa.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Kappa_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = ZEROMATRIX(3,3);


   return beta_Kappa;
}

/**
 * Calculates the 4-loop beta function of Kappa.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Kappa_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = ZEROMATRIX(3,3);


   return beta_Kappa;
}

/**
 * Calculates the 5-loop beta function of Kappa.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Kappa_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Kappa;

   beta_Kappa = ZEROMATRIX(3,3);


   return beta_Kappa;
}

} // namespace flexiblesusy
