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

// File generated at Fri 26 Jun 2015 18:59:15

#include "MRSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

/**
 * Calculates the one-loop beta function of mRd2.
 *
 * @return one-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRd2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mRd2;

   beta_mRd2 = Re(oneOver16PiSqr*(0.7745966692414834*g1*Tr11 + 2*(mHd2 +
      mRd2 + mS2)*AbsSqr(LamSD) + 3*(mHd2 + mRd2 + mT2)*AbsSqr(LamTD)));


   return beta_mRd2;
}

/**
 * Calculates the two-loop beta function of mRd2.
 *
 * @return two-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRd2_two_loop(const Soft_traces& soft_traces) const
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


   double beta_mRd2;

   beta_mRd2 = Re(twoLoop*(6*Power(g2,4)*Tr22 + 3.0983866769659336*g1*
      Tr31 - 2*AbsSqr(LamSD)*(3*tracemd2YdAdjYd + traceme2YeAdjYe +
      traceml2AdjYeYe + 3*tracemq2AdjYdYd + 6*mHd2*traceYdAdjYd + 3*mRd2*
      traceYdAdjYd + 3*mS2*traceYdAdjYd + 2*mHd2*traceYeAdjYe + mRd2*
      traceYeAdjYe + mS2*traceYeAdjYe + 2*(mHd2 + mHu2 + mRd2 + mRu2 + 2*mS2)*
      AbsSqr(LamSU) + 3*(2*mHd2 + 2*mRd2 + mS2 + mT2)*AbsSqr(LamTD)) + 1.2*
      Tr2U111*Sqr(g1) + 3*AbsSqr(LamTD)*(-3*tracemd2YdAdjYd - traceme2YeAdjYe -
      traceml2AdjYeYe - 3*tracemq2AdjYdYd - 3*(2*mHd2 + mRd2 + mT2)*
      traceYdAdjYd - 2*mHd2*traceYeAdjYe - mRd2*traceYeAdjYe - mT2*traceYeAdjYe
      - (mHd2 + mHu2 + mRd2 + mRu2 + 2*mT2)*AbsSqr(LamTU) + 4*mHd2*Sqr(g2) + 4
      *mRd2*Sqr(g2) + 4*mT2*Sqr(g2)) - 12*(mHd2 + mRd2 + mS2)*Sqr(LamSD)*Sqr(
      Conj(LamSD)) - 15*(mHd2 + mRd2 + mT2)*Sqr(LamTD)*Sqr(Conj(LamTD))));


   return beta_mRd2;
}

/**
 * Calculates the three-loop beta function of mRd2.
 *
 * @return three-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRd2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRd2;

   beta_mRd2 = 0;


   return beta_mRd2;
}

} // namespace flexiblesusy
