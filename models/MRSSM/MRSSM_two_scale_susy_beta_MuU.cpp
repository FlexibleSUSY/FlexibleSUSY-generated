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

// File generated at Sun 31 May 2015 12:24:50

#include "MRSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of MuU.
 *
 * @return one-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuU_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MuU;

   beta_MuU = Re(oneOver16PiSqr*(2*MuU*AbsSqr(LamSU) - 0.6*MuU*(-5*
      traceYuAdjYu - 5*AbsSqr(LamTU) + Sqr(g1) + 5*Sqr(g2))));


   return beta_MuU;
}

/**
 * Calculates the two-loop beta function of MuU.
 *
 * @return two-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuU_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_MuU;

   beta_MuU = Re(0.1*MuU*twoLoop*(45*Power(g1,4) + 165*Power(g2,4) - 30*
      traceYdAdjYuYuAdjYd - 90*traceYuAdjYuYuAdjYu - 40*AbsSqr(LamSD)*AbsSqr(
      LamSU) - 45*traceYuAdjYu*AbsSqr(LamTU) - 30*AbsSqr(LamTD)*AbsSqr(LamTU) -
      30*AbsSqr(LamSU)*(traceYuAdjYu + 2*AbsSqr(LamTU)) + 8*traceYuAdjYu*Sqr(
      g1) + 120*AbsSqr(LamTU)*Sqr(g2) + 18*Sqr(g1)*Sqr(g2) + 160*traceYuAdjYu*
      Sqr(g3) - 60*Sqr(LamSU)*Sqr(Conj(LamSU)) - 75*Sqr(LamTU)*Sqr(Conj(LamTU))
      ));


   return beta_MuU;
}

} // namespace flexiblesusy
