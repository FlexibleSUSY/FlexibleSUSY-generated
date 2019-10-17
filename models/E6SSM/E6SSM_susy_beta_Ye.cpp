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

// File generated at Wed 16 Oct 2019 21:52:05

#include "E6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Ye.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Ye_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(-0.1*Ye*(-30*traceYdAdjYd - 10*traceYeAdjYe - 10*
      AbsSqr(Lambdax) + 18*Sqr(g1) + 30*Sqr(g2) + 7*Sqr(gN)) + 3*(Ye*Ye.adjoint
      ()*Ye))).real();


   return beta_Ye;
}

/**
 * Calculates the 2-loop beta function of Ye.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Ye_2_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(0.025*Ye*(-360*traceYdAdjYdYdAdjYd - 120*
      traceYdAdjYuYuAdjYd - 120*traceYeAdjYeYeAdjYe - 120*traceKappaAdjKappa*
      AbsSqr(Lambdax) - 80*traceLambda12AdjLambda12*AbsSqr(Lambdax) - 120*
      traceYuAdjYu*AbsSqr(Lambdax) + 756*Quad(g1) + 660*Quad(g2) + 273*Quad(gN)
      - 16*traceYdAdjYd*Sqr(g1) + 48*traceYeAdjYe*Sqr(g1) + 72*Sqr(g1)*Sqr(g2)
      + 640*traceYdAdjYd*Sqr(g3) - 24*traceYdAdjYd*Sqr(gN) - 8*traceYeAdjYe*Sqr
      (gN) + 40*AbsSqr(Lambdax)*Sqr(gN) + 6*Sqr(g1)*Sqr(gN) + 78*Sqr(g2)*Sqr(gN
      ) - 120*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + 1.5*(-6*traceYdAdjYd - 2*
      traceYeAdjYe - 2*AbsSqr(Lambdax) + 4*Sqr(g2) + Sqr(gN))*(Ye*Ye.adjoint()*
      Ye) - 4*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();


   return beta_Ye;
}

/**
 * Calculates the 3-loop beta function of Ye.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Ye_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

/**
 * Calculates the 4-loop beta function of Ye.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Ye_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

/**
 * Calculates the 5-loop beta function of Ye.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Ye_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy
