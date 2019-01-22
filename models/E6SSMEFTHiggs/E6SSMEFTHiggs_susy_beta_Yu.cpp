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

// File generated at Tue 22 Jan 2019 14:41:49

#include "E6SSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Yu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(-0.03333333333333333*Yu*(-90*traceYuAdjYu - 30*
      AbsSqr(Lambdax) + 26*Sqr(g1) + 90*Sqr(g2) + 160*Sqr(g3) + 9*Sqr(gN)) + Yu
      *Yd.adjoint()*Yd + 3*(Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 2-loop beta function of Yu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Yu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0005555555555555556*Yu*(-5400*traceYdAdjYuYuAdjYd -
      16200*traceYuAdjYuYuAdjYu - 5400*traceKappaAdjKappa*AbsSqr(Lambdax) -
      3600*traceLambda12AdjLambda12*AbsSqr(Lambdax) - 5400*traceYdAdjYd*AbsSqr(
      Lambdax) - 1800*traceYeAdjYe*AbsSqr(Lambdax) + 15652*Quad(g1) + 29700*
      Quad(g2) + 25600*Quad(g3) + 5157*Quad(gN) + 1440*traceYuAdjYu*Sqr(g1) +
      1800*Sqr(g1)*Sqr(g2) + 28800*traceYuAdjYu*Sqr(g3) + 5440*Sqr(g1)*Sqr(g3)
      + 14400*Sqr(g2)*Sqr(g3) - 540*traceYuAdjYu*Sqr(gN) + 2700*AbsSqr(Lambdax)
      *Sqr(gN) + 966*Sqr(g1)*Sqr(gN) + 1350*Sqr(g2)*Sqr(gN) + 960*Sqr(g3)*Sqr(
      gN) - 5400*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + 0.2*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 5*AbsSqr(Lambdax) + 2*Sqr(g1) + 3*Sqr(gN))*(Yu*Yd.adjoint(
      )*Yd) + 0.2*(-45*traceYuAdjYu - 15*AbsSqr(Lambdax) + 2*Sqr(g1) + 30*Sqr(
      g2) + 3*Sqr(gN))*(Yu*Yu.adjoint()*Yu) - 2*(Yu*Yd.adjoint()*Yd*Yd.adjoint(
      )*Yd) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 4*(Yu*Yu.adjoint()*Yu*Yu
      .adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 3-loop beta function of Yu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Yu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

/**
 * Calculates the 4-loop beta function of Yu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Yu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

/**
 * Calculates the 5-loop beta function of Yu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_susy_parameters::calc_beta_Yu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
