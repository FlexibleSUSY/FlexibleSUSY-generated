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

// File generated at Wed 16 Oct 2019 21:40:22

#include "MRSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of MuD.
 *
 * @return 1-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuD_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_MuD;

   beta_MuD = Re(-0.2*MuD*oneOver16PiSqr*(-15*traceYdAdjYd - 5*traceYeAdjYe -
      10*AbsSqr(LamSD) - 15*AbsSqr(LamTD) + 3*Sqr(g1) + 15*Sqr(g2)));


   return beta_MuD;
}

/**
 * Calculates the 2-loop beta function of MuD.
 *
 * @return 2-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuD_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_MuD;

   beta_MuD = Re(0.1*MuD*twoLoop*(-90*traceYdAdjYdYdAdjYd - 30*
      traceYdAdjYuYuAdjYd - 30*traceYeAdjYeYeAdjYe - 30*traceYdAdjYd*AbsSqr(
      LamSD) - 10*traceYeAdjYe*AbsSqr(LamSD) - 40*AbsSqr(LamSD)*AbsSqr(LamSU) -
      45*traceYdAdjYd*AbsSqr(LamTD) - 15*traceYeAdjYe*AbsSqr(LamTD) - 60*AbsSqr
      (LamSD)*AbsSqr(LamTD) - 30*AbsSqr(LamTD)*AbsSqr(LamTU) + 45*Quad(g1) +
      165*Quad(g2) - 4*traceYdAdjYd*Sqr(g1) + 12*traceYeAdjYe*Sqr(g1) + 120*
      AbsSqr(LamTD)*Sqr(g2) + 18*Sqr(g1)*Sqr(g2) + 160*traceYdAdjYd*Sqr(g3) -
      60*Sqr(LamSD)*Sqr(Conj(LamSD)) - 75*Sqr(LamTD)*Sqr(Conj(LamTD))));


   return beta_MuD;
}

/**
 * Calculates the 3-loop beta function of MuD.
 *
 * @return 3-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuD_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuD;

   beta_MuD = 0;


   return beta_MuD;
}

/**
 * Calculates the 4-loop beta function of MuD.
 *
 * @return 4-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuD_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuD;

   beta_MuD = 0;


   return beta_MuD;
}

/**
 * Calculates the 5-loop beta function of MuD.
 *
 * @return 5-loop beta function
 */
double MRSSM_susy_parameters::calc_beta_MuD_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MuD;

   beta_MuD = 0;


   return beta_MuD;
}

} // namespace flexiblesusy
