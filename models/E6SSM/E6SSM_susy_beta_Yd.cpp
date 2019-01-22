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

// File generated at Tue 22 Jan 2019 17:01:19

#include "E6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Yd_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(-0.03333333333333333*Yd*(-90*traceYdAdjYd - 30*
      traceYeAdjYe - 30*AbsSqr(Lambdax) + 14*Sqr(g1) + 90*Sqr(g2) + 160*Sqr(g3)
      + 21*Sqr(gN)) + 3*(Yd*Yd.adjoint()*Yd) + Yd*Yu.adjoint()*Yu)).real();


   return beta_Yd;
}

/**
 * Calculates the 2-loop beta function of Yd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Yd_2_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(0.002777777777777778*Yd*(-3240*traceYdAdjYdYdAdjYd -
      1080*traceYdAdjYuYuAdjYd - 1080*traceYeAdjYeYeAdjYe - 1080*
      traceKappaAdjKappa*AbsSqr(Lambdax) - 720*traceLambda12AdjLambda12*AbsSqr(
      Lambdax) - 1080*traceYuAdjYu*AbsSqr(Lambdax) + 1652*Quad(g1) + 5940*Quad(
      g2) + 5120*Quad(g3) + 2457*Quad(gN) - 144*traceYdAdjYd*Sqr(g1) + 432*
      traceYeAdjYe*Sqr(g1) + 360*Sqr(g1)*Sqr(g2) + 5760*traceYdAdjYd*Sqr(g3) +
      320*Sqr(g1)*Sqr(g3) + 2880*Sqr(g2)*Sqr(g3) - 216*traceYdAdjYd*Sqr(gN) -
      72*traceYeAdjYe*Sqr(gN) + 360*AbsSqr(Lambdax)*Sqr(gN) - 84*Sqr(g1)*Sqr(gN
      ) + 540*Sqr(g2)*Sqr(gN) + 480*Sqr(g3)*Sqr(gN) - 1080*Sqr(Conj(Lambdax))*
      Sqr(Lambdax)) + 0.2*(-45*traceYdAdjYd - 15*traceYeAdjYe - 15*AbsSqr(
      Lambdax) + 4*Sqr(g1) + 30*Sqr(g2) + 6*Sqr(gN))*(Yd*Yd.adjoint()*Yd) + 0.2
      *(-15*traceYuAdjYu - 5*AbsSqr(Lambdax) + 4*Sqr(g1) + Sqr(gN))*(Yd*Yu.
      adjoint()*Yu) - 4*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint
      ()*Yu*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the 3-loop beta function of Yd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Yd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

/**
 * Calculates the 4-loop beta function of Yd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Yd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

/**
 * Calculates the 5-loop beta function of Yd.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_susy_parameters::calc_beta_Yd_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy
