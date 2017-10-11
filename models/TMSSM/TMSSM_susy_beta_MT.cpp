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

// File generated at Tue 10 Oct 2017 21:34:39

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

   beta_MT = Re(oneOver16PiSqr*(2*MT*AbsSqr(Lambdax) - 8*MT*Sqr(g2)));


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

   beta_MT = Re(0.4*MT*twoLoop*(140*Quad(g2) + AbsSqr(Lambdax)*(3*Sqr(g1)
      - 5*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + Sqr(g2))) - 15*Sqr
      (Conj(Lambdax))*Sqr(Lambdax)));


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

} // namespace flexiblesusy
