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

// File generated at Sun 4 Aug 2019 19:24:34

#include "MRSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of LamTU.
 *
 * @return 1-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamTU_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_LamTU;

   beta_LamTU = Re(0.2*LamTU*oneOver16PiSqr*(15*traceYuAdjYu + 10*AbsSqr(LamSU)
      + 5*AbsSqr(LamTD) + 20*AbsSqr(LamTU) - 3*Sqr(g1) - 35*Sqr(g2)));


   return beta_LamTU;
}

/**
 * Calculates the 2-loop beta function of LamTU.
 *
 * @return 2-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamTU_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_LamTU;

   beta_LamTU = Re(-0.1*LamTU*twoLoop*(30*traceYdAdjYuYuAdjYd + 90*
      traceYuAdjYuYuAdjYu + 30*traceYuAdjYu*AbsSqr(LamSU) + 40*AbsSqr(LamSD)*
      AbsSqr(LamSU) + 30*traceYdAdjYd*AbsSqr(LamTD) + 10*traceYeAdjYe*AbsSqr(
      LamTD) + 20*AbsSqr(LamSD)*AbsSqr(LamTD) + 75*traceYuAdjYu*AbsSqr(LamTU) +
      80*AbsSqr(LamSU)*AbsSqr(LamTU) + 30*AbsSqr(LamTD)*AbsSqr(LamTU) - 45*Quad
      (g1) - 485*Quad(g2) - 8*traceYuAdjYu*Sqr(g1) - 6*AbsSqr(LamTD)*Sqr(g1) -
      6*AbsSqr(LamTU)*Sqr(g1) + 10*AbsSqr(LamTD)*Sqr(g2) - 110*AbsSqr(LamTU)*
      Sqr(g2) - 18*Sqr(g1)*Sqr(g2) - 160*traceYuAdjYu*Sqr(g3) + 60*Sqr(LamSU)*
      Sqr(Conj(LamSU)) + 30*Sqr(LamTD)*Sqr(Conj(LamTD)) + 105*Sqr(LamTU)*Sqr(
      Conj(LamTU))));


   return beta_LamTU;
}

/**
 * Calculates the 3-loop beta function of LamTU.
 *
 * @return 3-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamTU_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamTU;

   beta_LamTU = 0;


   return beta_LamTU;
}

/**
 * Calculates the 4-loop beta function of LamTU.
 *
 * @return 4-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamTU_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamTU;

   beta_LamTU = 0;


   return beta_LamTU;
}

/**
 * Calculates the 5-loop beta function of LamTU.
 *
 * @return 5-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamTU_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamTU;

   beta_LamTU = 0;


   return beta_LamTU;
}

} // namespace flexiblesusy
