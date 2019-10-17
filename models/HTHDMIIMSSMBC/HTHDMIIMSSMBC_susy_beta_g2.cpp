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

// File generated at Wed 16 Oct 2019 21:21:42

#include "HTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_g2_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = Re(-2.3333333333333335*oneOver16PiSqr*Cube(g2));


   return beta_g2;
}

/**
 * Calculates the 2-loop beta function of g2.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_g2_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2;

   beta_g2 = Re(0.16666666666666666*twoLoop*Cube(g2)*(-9*traceYdAdjYd - 3*
      traceYeAdjYe - 9*traceYuAdjYu + 9*Sqr(g1) + 97*Sqr(g2) + 72*Sqr(g3)));


   return beta_g2;
}

/**
 * Calculates the 3-loop beta function of g2.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_g2_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return beta_g2;
}

/**
 * Calculates the 4-loop beta function of g2.
 *
 * @return 4-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_g2_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return beta_g2;
}

/**
 * Calculates the 5-loop beta function of g2.
 *
 * @return 5-loop beta function
 */
double HTHDMIIMSSMBC_susy_parameters::calc_beta_g2_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return beta_g2;
}

} // namespace flexiblesusy
