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

// File generated at Fri 10 Apr 2020 19:46:22

#include "HSSUSY_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g1.
 *
 * @return 1-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g1_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g1;

   beta_g1 = Re(4.1*oneOver16PiSqr*Cube(g1));


   return beta_g1;
}

/**
 * Calculates the 2-loop beta function of g1.
 *
 * @return 2-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g1_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g1;

   beta_g1 = Re(0.02*twoLoop*Cube(g1)*(-25*traceYdAdjYd - 75*traceYeAdjYe - 85*
      traceYuAdjYu + 199*Sqr(g1) + 135*Sqr(g2) + 440*Sqr(g3)));


   return beta_g1;
}

/**
 * Calculates the 3-loop beta function of g1.
 *
 * @return 3-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1;

   beta_g1 = Re(-0.000041666666666666665*threeLoop*Cube(g1)*(388613*Quad(g1) -
      295875*Quad(g2) - 1425600*Quad(g3) - 283500*Quad(Yu(2,2)) - 6480*Lambdax*
      Sqr(g1) - 10800*Lambdax*Sqr(g2) - 18450*Sqr(g1)*Sqr(g2) + 43840*Sqr(g1)*
      Sqr(g3) + 14400*Sqr(g2)*Sqr(g3) + 10800*Sqr(Lambdax) + 84810*Sqr(g1)*Sqr(
      Yu(2,2)) + 353250*Sqr(g2)*Sqr(Yu(2,2)) + 139200*Sqr(g3)*Sqr(Yu(2,2))));


   return beta_g1;
}

/**
 * Calculates the 4-loop beta function of g1.
 *
 * @return 4-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g1_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1;

   beta_g1 = 0;


   return beta_g1;
}

/**
 * Calculates the 5-loop beta function of g1.
 *
 * @return 5-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g1_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1;

   beta_g1 = 0;


   return beta_g1;
}

} // namespace flexiblesusy
