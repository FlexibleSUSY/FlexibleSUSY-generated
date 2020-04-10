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

// File generated at Fri 10 Apr 2020 19:53:54

#include "TMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of MT.
 *
 * @return 1-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_MT_1_loop(const Susy_traces& susy_traces) const
{


   double beta_MT;

   beta_MT = Re(-2*MT*oneOver16PiSqr*(-AbsSqr(Lambdax) + 4*Sqr(g2)));


   return beta_MT;
}

/**
 * Calculates the 2-loop beta function of MT.
 *
 * @return 2-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_MT_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MT;

   beta_MT = Re(0.4*MT*twoLoop*(-15*traceYdAdjYd*AbsSqr(Lambdax) - 5*
      traceYeAdjYe*AbsSqr(Lambdax) - 15*traceYuAdjYu*AbsSqr(Lambdax) + 140*Quad
      (g2) + 3*AbsSqr(Lambdax)*Sqr(g1) - 5*AbsSqr(Lambdax)*Sqr(g2) - 15*Sqr(
      Conj(Lambdax))*Sqr(Lambdax)));


   return beta_MT;
}

/**
 * Calculates the 3-loop beta function of MT.
 *
 * @return 3-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_MT_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MT;

   beta_MT = 0;


   return beta_MT;
}

/**
 * Calculates the 4-loop beta function of MT.
 *
 * @return 4-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_MT_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MT;

   beta_MT = 0;


   return beta_MT;
}

/**
 * Calculates the 5-loop beta function of MT.
 *
 * @return 5-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_MT_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MT;

   beta_MT = 0;


   return beta_MT;
}

} // namespace flexiblesusy
