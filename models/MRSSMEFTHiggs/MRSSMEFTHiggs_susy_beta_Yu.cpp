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

// File generated at Sun 26 Aug 2018 14:25:04

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

   beta_Yu = (oneOver16PiSqr*(Yu*(3*traceYuAdjYu + AbsSqr(LamSU) + 1.5*AbsSqr(
      LamTU) - 0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(
      g3)) + Yu*Yd.adjoint()*Yd + 3*(Yu*Yu.adjoint()*Yu))).real();


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

   beta_Yu = (twoLoop*(Yu*(-3*traceYdAdjYuYuAdjYd - 9*traceYuAdjYuYuAdjYu - 2*
      AbsSqr(LamSD)*AbsSqr(LamSU) - 3*AbsSqr(LamSU)*AbsSqr(LamTU) - 1.5*AbsSqr(
      LamTD)*AbsSqr(LamTU) + 6.615555555555556*Quad(g1) + 16.5*Quad(g2) +
      14.222222222222221*Quad(g3) + 0.8*traceYuAdjYu*Sqr(g1) + 6*AbsSqr(LamTU)*
      Sqr(g2) + Sqr(g1)*Sqr(g2) + 16*traceYuAdjYu*Sqr(g3) + 3.022222222222222*
      Sqr(g1)*Sqr(g3) + 8*Sqr(g2)*Sqr(g3) - 3*Sqr(LamSU)*Sqr(Conj(LamSU)) -
      3.75*Sqr(LamTU)*Sqr(Conj(LamTU))) + (-3*traceYdAdjYd - traceYeAdjYe -
      AbsSqr(LamSD) - 1.5*AbsSqr(LamTD) + 0.4*Sqr(g1))*(Yu*Yd.adjoint()*Yd) + (
      -9*traceYuAdjYu - 3*AbsSqr(LamSU) - 4.5*AbsSqr(LamTU) + 0.4*Sqr(g1) + 6*
      Sqr(g2))*(Yu*Yu.adjoint()*Yu) - 2*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) -
      2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 4*(Yu*Yu.adjoint()*Yu*Yu.adjoint
      ()*Yu))).real();


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

} // namespace flexiblesusy
