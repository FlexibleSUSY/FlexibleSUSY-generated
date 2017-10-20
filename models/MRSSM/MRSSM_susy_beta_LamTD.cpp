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

// File generated at Fri 20 Oct 2017 08:45:08

#include "MRSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of LamTD.
 *
 * @return 1-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamTD_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_LamTD;

   beta_LamTD = Re(oneOver16PiSqr*(3*LamTD*traceYdAdjYd + LamTD*
      traceYeAdjYe + 2*LamTD*AbsSqr(LamSD) + LamTD*AbsSqr(LamTU) - 0.6*LamTD*
      Sqr(g1) - 7*LamTD*Sqr(g2) + 4*Conj(LamTD)*Sqr(LamTD)));


   return beta_LamTD;
}

/**
 * Calculates the 2-loop beta function of LamTD.
 *
 * @return 2-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamTD_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_LamTD;

   beta_LamTD = Re(-0.1*LamTD*twoLoop*(90*traceYdAdjYdYdAdjYd + 30*
      traceYdAdjYuYuAdjYd + 30*traceYeAdjYeYeAdjYe + 10*AbsSqr(LamSD)*(3*
      traceYdAdjYd + traceYeAdjYe + 4*AbsSqr(LamSU) + 8*AbsSqr(LamTD)) + 30*
      traceYuAdjYu*AbsSqr(LamTU) + 20*AbsSqr(LamSU)*AbsSqr(LamTU) - 45*Quad(g1)
      - 485*Quad(g2) + 4*traceYdAdjYd*Sqr(g1) - 12*traceYeAdjYe*Sqr(g1) - 6*
      AbsSqr(LamTU)*Sqr(g1) + AbsSqr(LamTD)*(75*traceYdAdjYd + 25*traceYeAdjYe
      + 30*AbsSqr(LamTU) - 6*Sqr(g1) - 110*Sqr(g2)) + 10*AbsSqr(LamTU)*Sqr(g2)
      - 18*Sqr(g1)*Sqr(g2) - 160*traceYdAdjYd*Sqr(g3) + 60*Sqr(LamSD)*Sqr(Conj(
      LamSD)) + 105*Sqr(LamTD)*Sqr(Conj(LamTD)) + 30*Sqr(LamTU)*Sqr(Conj(LamTU)
      )));


   return beta_LamTD;
}

/**
 * Calculates the 3-loop beta function of LamTD.
 *
 * @return 3-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamTD_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamTD;

   beta_LamTD = 0;


   return beta_LamTD;
}

} // namespace flexiblesusy
