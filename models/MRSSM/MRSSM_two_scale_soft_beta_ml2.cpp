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

// File generated at Sun 18 Oct 2015 11:51:11

#include "MRSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of ml2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_ml2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (oneOver16PiSqr*(2*mHd2*(Ye.adjoint()*Ye) + ml2*Ye.adjoint(
      )*Ye + 2*(Ye.adjoint()*me2*Ye) + Ye.adjoint()*Ye*ml2 - 0.7745966692414834
      *g1*Tr11*UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the two-loop beta function of ml2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_ml2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (twoLoop*(-0.2*(10*(2*mHd2 + mRd2 + mS2)*AbsSqr(LamSD) + 15
      *(2*mHd2 + mRd2 + mT2)*AbsSqr(LamTD) + 2*(15*tracemd2YdAdjYd + 5*
      traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*tracemq2AdjYdYd + 30*mHd2*
      traceYdAdjYd + 10*mHd2*traceYeAdjYe - 6*mHd2*Sqr(g1)))*(Ye.adjoint()*Ye)
      - 3*traceYdAdjYd*(ml2*Ye.adjoint()*Ye) - traceYeAdjYe*(ml2*Ye.adjoint()*
      Ye) - AbsSqr(LamSD)*(ml2*Ye.adjoint()*Ye) - 1.5*AbsSqr(LamTD)*(ml2*
      Ye.adjoint()*Ye) + 1.2*Sqr(g1)*(ml2*Ye.adjoint()*Ye) - 6*traceYdAdjYd*(
      Ye.adjoint()*me2*Ye) - 2*traceYeAdjYe*(Ye.adjoint()*me2*Ye) - 2*AbsSqr(
      LamSD)*(Ye.adjoint()*me2*Ye) - 3*AbsSqr(LamTD)*(Ye.adjoint()*me2*Ye) +
      2.4*Sqr(g1)*(Ye.adjoint()*me2*Ye) - 3*traceYdAdjYd*(Ye.adjoint()*Ye*ml2)
      - traceYeAdjYe*(Ye.adjoint()*Ye*ml2) - AbsSqr(LamSD)*(Ye.adjoint()*Ye*ml2
      ) - 1.5*AbsSqr(LamTD)*(Ye.adjoint()*Ye*ml2) + 1.2*Sqr(g1)*(Ye.adjoint()*
      Ye*ml2) - 8*mHd2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 2*(ml2*Ye.adjoint()*
      Ye*Ye.adjoint()*Ye) - 4*(Ye.adjoint()*me2*Ye*Ye.adjoint()*Ye) - 4*(
      Ye.adjoint()*Ye*ml2*Ye.adjoint()*Ye) - 4*(Ye.adjoint()*Ye*Ye.adjoint()*
      me2*Ye) - 2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye*ml2) + (6*Power(g2,4)*Tr22 +
      0.4*g1*(3*g1*Tr2U111 - 7.745966692414834*Tr31))*UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the three-loop beta function of ml2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_ml2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

} // namespace flexiblesusy
