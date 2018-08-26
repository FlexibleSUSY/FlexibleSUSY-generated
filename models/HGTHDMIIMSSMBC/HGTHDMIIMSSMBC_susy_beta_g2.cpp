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

// File generated at Sun 26 Aug 2018 14:06:05

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = Re(-(oneOver16PiSqr*Cube(g2)));


   return beta_g2;
}

/**
 * Calculates the 2-loop beta function of g2.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2;

   beta_g2 = Re(0.25*twoLoop*Cube(g2)*(-6*traceYdAdjYd - 2*traceYeAdjYe - 6*
      traceYuAdjYu + 6*Sqr(g1) - 11*Sqr(g1d) - Sqr(g1dp) + 150*Sqr(g2) - 11*Sqr
      (g2u) - Sqr(g2up) + 48*Sqr(g3)));


   return beta_g2;
}

/**
 * Calculates the 3-loop beta function of g2.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2_3_loop(const Susy_traces& susy_traces) const
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
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return beta_g2;
}

} // namespace flexiblesusy
