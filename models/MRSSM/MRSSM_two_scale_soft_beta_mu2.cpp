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

// File generated at Mon 23 Feb 2015 12:39:12

#include "MRSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mu2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_mu2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = oneOver16PiSqr*(2*(2*mHu2*(Yu*Yu.adjoint()) + mu2*Yu*
      Yu.adjoint() + 2*(Yu*mq2*Yu.adjoint()) + Yu*Yu.adjoint()*mu2) -
      1.0327955589886444*g1*Tr11*UNITMATRIX(3));


   return beta_mu2;
}

/**
 * Calculates the two-loop beta function of mu2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_mu2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = twoLoop*(-0.4*(30*tracemq2AdjYuYu + 30*tracemu2YuAdjYu + 60
      *mHu2*traceYuAdjYu + 10*(2*mHu2 + mRu2 + mS2)*AbsSqr(LamSU) + 15*(2*mHu2
      + mRu2 + mT2)*AbsSqr(LamTU) + 2*mHu2*Sqr(g1) - 30*mHu2*Sqr(g2))*(Yu*
      Yu.adjoint()) - 6*traceYuAdjYu*(mu2*Yu*Yu.adjoint()) - 2*AbsSqr(LamSU)*(
      mu2*Yu*Yu.adjoint()) - 3*AbsSqr(LamTU)*(mu2*Yu*Yu.adjoint()) - 0.4*Sqr(g1
      )*(mu2*Yu*Yu.adjoint()) + 6*Sqr(g2)*(mu2*Yu*Yu.adjoint()) - 12*
      traceYuAdjYu*(Yu*mq2*Yu.adjoint()) - 4*AbsSqr(LamSU)*(Yu*mq2*Yu.adjoint()
      ) - 6*AbsSqr(LamTU)*(Yu*mq2*Yu.adjoint()) - 0.8*Sqr(g1)*(Yu*mq2*
      Yu.adjoint()) + 12*Sqr(g2)*(Yu*mq2*Yu.adjoint()) - 6*traceYuAdjYu*(Yu*
      Yu.adjoint()*mu2) - 2*AbsSqr(LamSU)*(Yu*Yu.adjoint()*mu2) - 3*AbsSqr(
      LamTU)*(Yu*Yu.adjoint()*mu2) - 0.4*Sqr(g1)*(Yu*Yu.adjoint()*mu2) + 6*Sqr(
      g2)*(Yu*Yu.adjoint()*mu2) - 4*mHd2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 4*
      mHu2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 8*mHu2*(Yu*Yu.adjoint()*Yu*
      Yu.adjoint()) - 2*(mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 2*(mu2*Yu*
      Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*mq2*Yd.adjoint()*Yd*Yu.adjoint()) -
      4*(Yu*mq2*Yu.adjoint()*Yu*Yu.adjoint()) - 4*(Yu*Yd.adjoint()*md2*Yd*
      Yu.adjoint()) - 4*(Yu*Yd.adjoint()*Yd*mq2*Yu.adjoint()) - 2*(Yu*
      Yd.adjoint()*Yd*Yu.adjoint()*mu2) - 4*(Yu*Yu.adjoint()*mu2*Yu*Yu.adjoint(
      )) - 4*(Yu*Yu.adjoint()*Yu*mq2*Yu.adjoint()) - 2*(Yu*Yu.adjoint()*Yu*
      Yu.adjoint()*mu2) + 1.0666666666666667*(10*Power(g3,4)*Tr23 + g1*(2*g1*
      Tr2U111 - 3.872983346207417*Tr31))*UNITMATRIX(3));


   return beta_mu2;
}

} // namespace flexiblesusy
