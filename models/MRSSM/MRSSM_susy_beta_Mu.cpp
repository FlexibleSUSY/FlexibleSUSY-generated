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

// File generated at Sun 4 Aug 2019 19:24:37

#include "MRSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Mu.
 *
 * @return 1-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_Mu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(-0.1*oneOver16PiSqr*Mu*(-30*traceYdAdjYd - 10*traceYeAdjYe - 30
      *traceYuAdjYu - 10*AbsSqr(LamSD) - 10*AbsSqr(LamSU) - 15*AbsSqr(LamTD) -
      15*AbsSqr(LamTU) + 6*Sqr(g1) + 30*Sqr(g2)));


   return beta_Mu;
}

/**
 * Calculates the 2-loop beta function of Mu.
 *
 * @return 2-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_Mu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(0.05*twoLoop*Mu*(-180*traceYdAdjYdYdAdjYd - 120*
      traceYdAdjYuYuAdjYd - 60*traceYeAdjYeYeAdjYe - 180*traceYuAdjYuYuAdjYu -
      80*AbsSqr(LamSD)*AbsSqr(LamSU) - 60*AbsSqr(LamSD)*AbsSqr(LamTD) - 60*
      AbsSqr(LamSU)*AbsSqr(LamTU) - 60*AbsSqr(LamTD)*AbsSqr(LamTU) + 90*Quad(g1
      ) + 330*Quad(g2) - 8*traceYdAdjYd*Sqr(g1) + 24*traceYeAdjYe*Sqr(g1) + 16*
      traceYuAdjYu*Sqr(g1) + 120*AbsSqr(LamTD)*Sqr(g2) + 120*AbsSqr(LamTU)*Sqr(
      g2) + 36*Sqr(g1)*Sqr(g2) + 320*traceYdAdjYd*Sqr(g3) + 320*traceYuAdjYu*
      Sqr(g3) - 60*Sqr(LamSD)*Sqr(Conj(LamSD)) - 60*Sqr(LamSU)*Sqr(Conj(LamSU))
      - 75*Sqr(LamTD)*Sqr(Conj(LamTD)) - 75*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_Mu;
}

/**
 * Calculates the 3-loop beta function of Mu.
 *
 * @return 3-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_Mu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

/**
 * Calculates the 4-loop beta function of Mu.
 *
 * @return 4-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_Mu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

/**
 * Calculates the 5-loop beta function of Mu.
 *
 * @return 5-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_Mu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

} // namespace flexiblesusy
