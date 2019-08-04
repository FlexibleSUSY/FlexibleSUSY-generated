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

// File generated at Sun 4 Aug 2019 17:35:23

#include "MRSSMEFTHiggs_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of LamSD.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSD_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_LamSD;

   beta_LamSD = Re(0.2*LamSD*oneOver16PiSqr*(15*traceYdAdjYd + 5*traceYeAdjYe +
      20*AbsSqr(LamSD) + 10*AbsSqr(LamSU) + 15*AbsSqr(LamTD) - 3*Sqr(g1) - 15*
      Sqr(g2)));


   return beta_LamSD;
}

/**
 * Calculates the 2-loop beta function of LamSD.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSD_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_LamSD;

   beta_LamSD = Re(-0.1*LamSD*twoLoop*(90*traceYdAdjYdYdAdjYd + 30*
      traceYdAdjYuYuAdjYd + 30*traceYeAdjYeYeAdjYe + 90*traceYdAdjYd*AbsSqr(
      LamSD) + 30*traceYeAdjYe*AbsSqr(LamSD) + 60*traceYuAdjYu*AbsSqr(LamSU) +
      40*AbsSqr(LamSD)*AbsSqr(LamSU) + 45*traceYdAdjYd*AbsSqr(LamTD) + 15*
      traceYeAdjYe*AbsSqr(LamTD) + 120*AbsSqr(LamSD)*AbsSqr(LamTD) + 60*AbsSqr(
      LamSU)*AbsSqr(LamTU) + 30*AbsSqr(LamTD)*AbsSqr(LamTU) - 45*Quad(g1) - 165
      *Quad(g2) + 4*traceYdAdjYd*Sqr(g1) - 12*traceYeAdjYe*Sqr(g1) - 12*AbsSqr(
      LamSD)*Sqr(g1) - 12*AbsSqr(LamSU)*Sqr(g1) - 60*AbsSqr(LamSD)*Sqr(g2) - 60
      *AbsSqr(LamSU)*Sqr(g2) - 120*AbsSqr(LamTD)*Sqr(g2) - 18*Sqr(g1)*Sqr(g2) -
      160*traceYdAdjYd*Sqr(g3) + 100*Sqr(LamSD)*Sqr(Conj(LamSD)) + 40*Sqr(LamSU
      )*Sqr(Conj(LamSU)) + 75*Sqr(LamTD)*Sqr(Conj(LamTD))));


   return beta_LamSD;
}

/**
 * Calculates the 3-loop beta function of LamSD.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSD_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSD;

   beta_LamSD = 0;


   return beta_LamSD;
}

/**
 * Calculates the 4-loop beta function of LamSD.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSD_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSD;

   beta_LamSD = 0;


   return beta_LamSD;
}

/**
 * Calculates the 5-loop beta function of LamSD.
 *
 * @return 5-loop beta function
 */
double MRSSMEFTHiggs_susy_parameters::calc_beta_LamSD_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSD;

   beta_LamSD = 0;


   return beta_LamSD;
}

} // namespace flexiblesusy
