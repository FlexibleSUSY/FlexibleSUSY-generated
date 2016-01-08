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

// File generated at Fri 8 Jan 2016 15:09:43

#include "MRSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Mu.
 *
 * @return one-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_Mu_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(oneOver16PiSqr*(3*traceYdAdjYd*Mu + traceYeAdjYe*Mu + 3*
      traceYuAdjYu*Mu + AbsSqr(LamSD)*Mu + AbsSqr(LamSU)*Mu + 1.5*AbsSqr(LamTD)
      *Mu + 1.5*AbsSqr(LamTU)*Mu - 0.6*Mu*Sqr(g1) - 3*Mu*Sqr(g2)));


   return beta_Mu;
}

/**
 * Calculates the two-loop beta function of Mu.
 *
 * @return two-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_Mu_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(0.05*twoLoop*Mu*(90*Power(g1,4) + 330*Power(g2,4) - 180*
      traceYdAdjYdYdAdjYd - 120*traceYdAdjYuYuAdjYd - 60*traceYeAdjYeYeAdjYe -
      180*traceYuAdjYuYuAdjYu - 20*AbsSqr(LamSD)*(4*AbsSqr(LamSU) + 3*AbsSqr(
      LamTD)) - 60*AbsSqr(LamSU)*AbsSqr(LamTU) - 60*AbsSqr(LamTD)*AbsSqr(LamTU)
      - 8*traceYdAdjYd*Sqr(g1) + 24*traceYeAdjYe*Sqr(g1) + 16*traceYuAdjYu*Sqr
      (g1) + 120*AbsSqr(LamTD)*Sqr(g2) + 120*AbsSqr(LamTU)*Sqr(g2) + 36*Sqr(g1)
      *Sqr(g2) + 320*traceYdAdjYd*Sqr(g3) + 320*traceYuAdjYu*Sqr(g3) - 60*Sqr(
      LamSD)*Sqr(Conj(LamSD)) - 60*Sqr(LamSU)*Sqr(Conj(LamSU)) - 75*Sqr(LamTD)*
      Sqr(Conj(LamTD)) - 75*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_Mu;
}

/**
 * Calculates the three-loop beta function of Mu.
 *
 * @return three-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_Mu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

} // namespace flexiblesusy
