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

// File generated at Wed 16 Oct 2019 19:26:31

#include "MRSSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_susy_parameters::calc_beta_Yu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(-0.03333333333333333*Yu*(-90*traceYuAdjYu - 30*
      AbsSqr(LamSU) - 45*AbsSqr(LamTU) + 26*Sqr(g1) + 90*Sqr(g2) + 160*Sqr(g3))
      + Yu*Yd.adjoint()*Yd + 3*(Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 2-loop beta function of Yu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_susy_parameters::calc_beta_Yu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0011111111111111111*Yu*(-2700*traceYdAdjYuYuAdjYd -
      8100*traceYuAdjYuYuAdjYu - 1800*AbsSqr(LamSD)*AbsSqr(LamSU) - 2700*AbsSqr
      (LamSU)*AbsSqr(LamTU) - 1350*AbsSqr(LamTD)*AbsSqr(LamTU) + 5954*Quad(g1)
      + 14850*Quad(g2) + 12800*Quad(g3) + 720*traceYuAdjYu*Sqr(g1) + 5400*
      AbsSqr(LamTU)*Sqr(g2) + 900*Sqr(g1)*Sqr(g2) + 14400*traceYuAdjYu*Sqr(g3)
      + 2720*Sqr(g1)*Sqr(g3) + 7200*Sqr(g2)*Sqr(g3) - 2700*Sqr(LamSU)*Sqr(Conj(
      LamSU)) - 3375*Sqr(LamTU)*Sqr(Conj(LamTU))) + 0.1*(-30*traceYdAdjYd - 10*
      traceYeAdjYe - 10*AbsSqr(LamSD) - 15*AbsSqr(LamTD) + 4*Sqr(g1))*(Yu*Yd.
      adjoint()*Yd) + 0.1*(-90*traceYuAdjYu - 30*AbsSqr(LamSU) - 45*AbsSqr(
      LamTU) + 4*Sqr(g1) + 60*Sqr(g2))*(Yu*Yu.adjoint()*Yu) - 2*(Yu*Yd.adjoint(
      )*Yd*Yd.adjoint()*Yd) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 4*(Yu*Yu
      .adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 3-loop beta function of Yu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_susy_parameters::calc_beta_Yu_3_loop(const Susy_traces& susy_traces) const
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
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_susy_parameters::calc_beta_Yu_4_loop(const Susy_traces& susy_traces) const
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
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_susy_parameters::calc_beta_Yu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
