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

// File generated at Tue 22 Jan 2019 16:49:51

#include "TMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vT.
 *
 * @return 1-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_vT_1_loop(const Susy_traces& susy_traces) const
{


   double beta_vT;

   beta_vT = Re(oneOver16PiSqr*vT*(-AbsSqr(Lambdax) + 4*Sqr(g2)));


   return beta_vT;
}

/**
 * Calculates the 2-loop beta function of vT.
 *
 * @return 2-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_vT_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vT;

   beta_vT = Re(-0.06666666666666667*twoLoop*vT*(-45*traceYdAdjYd*AbsSqr(
      Lambdax) - 15*traceYeAdjYe*AbsSqr(Lambdax) - 45*traceYuAdjYu*AbsSqr(
      Lambdax) + 190*Quad(g2) + 9*AbsSqr(Lambdax)*Sqr(g1) + 45*AbsSqr(Lambdax)*
      Sqr(g2) - 45*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


   return beta_vT;
}

/**
 * Calculates the 3-loop beta function of vT.
 *
 * @return 3-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_vT_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vT;

   beta_vT = 0;


   return beta_vT;
}

/**
 * Calculates the 4-loop beta function of vT.
 *
 * @return 4-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_vT_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vT;

   beta_vT = 0;


   return beta_vT;
}

/**
 * Calculates the 5-loop beta function of vT.
 *
 * @return 5-loop beta function
 */
double TMSSM_susy_parameters::calc_beta_vT_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vT;

   beta_vT = 0;


   return beta_vT;
}

} // namespace flexiblesusy
