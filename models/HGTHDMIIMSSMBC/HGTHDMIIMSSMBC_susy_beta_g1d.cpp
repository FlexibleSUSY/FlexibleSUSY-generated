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

// File generated at Mon 5 Mar 2018 17:08:17

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g1d.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_g1d;

   beta_g1d = Re(0.05*g1d*oneOver16PiSqr*(-9*Sqr(g1) + 5*(12*traceYdAdjYd
      + 4*traceYeAdjYe + 11*Sqr(g1d) + 3*Sqr(g1dp) - 33*Sqr(g2) + 2*Sqr(g2u)))
      );


   return beta_g1d;
}

/**
 * Calculates the 2-loop beta function of g1d.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_g1d;

   beta_g1d = Re(0.00125*twoLoop*(516*g1d*Quad(g1) + 5*g1d*Sqr(g1)*(100*
      traceYdAdjYd + 300*traceYeAdjYe + 435*Sqr(g1d) + 63*Sqr(g1dp) + 72*Sqr(g2
      ) - 42*Sqr(g2u)) - 25*(112*Power5(g1d) + 16*g1dp*g2u*g2up*(2*Lambda4 + 3*
      Sqr(g2)) + Cube(g1d)*(320*Lambda1 + 180*traceYdAdjYd + 60*traceYeAdjYe +
      118*Sqr(g1dp) - 875*Sqr(g2) + 42*Sqr(g2u)) + g1d*(10*Quad(g1dp) + Sqr(
      g1dp)*(64*Lambda1 + 36*traceYdAdjYd + 12*traceYeAdjYe - 111*Sqr(g2) + 14*
      Sqr(g2up)) + 2*(508*Quad(g2) + 21*Quad(g2u) - Sqr(g2)*(30*(3*traceYdAdjYd
      + traceYeAdjYe) + 43*Sqr(g2u)) + Sqr(g2u)*(32*Lambda3 - 16*Lambda4 + 36*
      traceYuAdjYu + 7*Sqr(g2up)) - 4*(4*Lambda3*Lambda4 - 27*
      traceYdAdjYdYdAdjYd - 9*traceYdAdjYuYuAdjYd - 9*traceYeAdjYeYeAdjYe + 80*
      traceYdAdjYd*Sqr(g3) + 24*Sqr(Lambda1) + 4*Sqr(Lambda3) + 4*Sqr(Lambda4)
      + 6*Sqr(Lambda5) + 18*Sqr(Lambda6) + 6*Sqr(Lambda7)))))));


   return beta_g1d;
}

/**
 * Calculates the 3-loop beta function of g1d.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1d;

   beta_g1d = 0;


   return beta_g1d;
}

/**
 * Calculates the 4-loop beta function of g1d.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g1d_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1d;

   beta_g1d = 0;


   return beta_g1d;
}

} // namespace flexiblesusy
