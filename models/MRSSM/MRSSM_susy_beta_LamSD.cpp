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

// File generated at Mon 5 Mar 2018 17:47:50

#include "MRSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of LamSD.
 *
 * @return 1-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamSD_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_LamSD;

   beta_LamSD = Re(oneOver16PiSqr*(3*LamSD*traceYdAdjYd + LamSD*
      traceYeAdjYe + 2*LamSD*AbsSqr(LamSU) + 3*LamSD*AbsSqr(LamTD) - 0.6*LamSD*
      Sqr(g1) - 3*LamSD*Sqr(g2) + 4*Conj(LamSD)*Sqr(LamSD)));


   return beta_LamSD;
}

/**
 * Calculates the 2-loop beta function of LamSD.
 *
 * @return 2-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamSD_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_LamSD;

   beta_LamSD = Re(-0.1*LamSD*twoLoop*(90*traceYdAdjYdYdAdjYd + 30*
      traceYdAdjYuYuAdjYd + 30*traceYeAdjYeYeAdjYe + 45*traceYdAdjYd*AbsSqr(
      LamTD) + 15*traceYeAdjYe*AbsSqr(LamTD) + 30*AbsSqr(LamTD)*AbsSqr(LamTU) -
      45*Quad(g1) - 165*Quad(g2) + 4*traceYdAdjYd*Sqr(g1) - 12*traceYeAdjYe*
      Sqr(g1) + 2*AbsSqr(LamSD)*(45*traceYdAdjYd + 15*traceYeAdjYe + 20*AbsSqr(
      LamSU) + 60*AbsSqr(LamTD) - 6*Sqr(g1) - 30*Sqr(g2)) - 120*AbsSqr(LamTD)*
      Sqr(g2) - 18*Sqr(g1)*Sqr(g2) - 12*AbsSqr(LamSU)*(-5*traceYuAdjYu - 5*
      AbsSqr(LamTU) + Sqr(g1) + 5*Sqr(g2)) - 160*traceYdAdjYd*Sqr(g3) + 100*Sqr
      (LamSD)*Sqr(Conj(LamSD)) + 40*Sqr(LamSU)*Sqr(Conj(LamSU)) + 75*Sqr(LamTD)
      *Sqr(Conj(LamTD))));


   return beta_LamSD;
}

/**
 * Calculates the 3-loop beta function of LamSD.
 *
 * @return 3-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_LamSD_3_loop(const Susy_traces& susy_traces) const
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
double MRSSM_susy_parameters::calc_beta_LamSD_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_LamSD;

   beta_LamSD = 0;


   return beta_LamSD;
}

} // namespace flexiblesusy
