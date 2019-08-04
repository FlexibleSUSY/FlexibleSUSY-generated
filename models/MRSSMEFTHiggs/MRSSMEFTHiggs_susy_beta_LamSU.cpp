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

// File generated at Sun 4 Aug 2019 17:35:24

#include "MRSSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of LamSU.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSU_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_LamSU;

   beta_LamSU = Re(0.2*LamSU*oneOver16PiSqr*(15*traceYuAdjYu + 10*AbsSqr(LamSD)
      + 20*AbsSqr(LamSU) + 15*AbsSqr(LamTU) - 3*Sqr(g1) - 15*Sqr(g2)));


   return beta_LamSU;
}

/**
 * Calculates the 2-loop beta function of LamSU.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSU_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_LamSU;

   beta_LamSU = Re(-0.1*LamSU*twoLoop*(30*traceYdAdjYuYuAdjYd + 90*
      traceYuAdjYuYuAdjYu + 60*traceYdAdjYd*AbsSqr(LamSD) + 20*traceYeAdjYe*
      AbsSqr(LamSD) + 90*traceYuAdjYu*AbsSqr(LamSU) + 40*AbsSqr(LamSD)*AbsSqr(
      LamSU) + 60*AbsSqr(LamSD)*AbsSqr(LamTD) + 45*traceYuAdjYu*AbsSqr(LamTU) +
      120*AbsSqr(LamSU)*AbsSqr(LamTU) + 30*AbsSqr(LamTD)*AbsSqr(LamTU) - 45*
      Quad(g1) - 165*Quad(g2) - 8*traceYuAdjYu*Sqr(g1) - 12*AbsSqr(LamSD)*Sqr(
      g1) - 12*AbsSqr(LamSU)*Sqr(g1) - 60*AbsSqr(LamSD)*Sqr(g2) - 60*AbsSqr(
      LamSU)*Sqr(g2) - 120*AbsSqr(LamTU)*Sqr(g2) - 18*Sqr(g1)*Sqr(g2) - 160*
      traceYuAdjYu*Sqr(g3) + 40*Sqr(LamSD)*Sqr(Conj(LamSD)) + 100*Sqr(LamSU)*
      Sqr(Conj(LamSU)) + 75*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_LamSU;
}

/**
 * Calculates the 3-loop beta function of LamSU.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSU_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSU;

   beta_LamSU = 0;


   return beta_LamSU;
}

/**
 * Calculates the 4-loop beta function of LamSU.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSU_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSU;

   beta_LamSU = 0;


   return beta_LamSU;
}

/**
 * Calculates the 5-loop beta function of LamSU.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSU_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSU;

   beta_LamSU = 0;


   return beta_LamSU;
}

} // namespace flexiblesusy
