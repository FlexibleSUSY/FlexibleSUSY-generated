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

// File generated at Fri 10 Apr 2020 19:47:19

#include "SM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g3.
 *
 * @return 1-loop beta function
 */
double SM_susy_parameters::calc_beta_g3_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g3;

   beta_g3 = Re(-7*oneOver16PiSqr*Cube(g3));


   return beta_g3;
}

/**
 * Calculates the 2-loop beta function of g3.
 *
 * @return 2-loop beta function
 */
double SM_susy_parameters::calc_beta_g3_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g3;

   beta_g3 = Re(-0.1*twoLoop*Cube(g3)*(20*traceYdAdjYd + 20*traceYuAdjYu - 11*
      Sqr(g1) - 45*Sqr(g2) + 260*Sqr(g3)));


   return beta_g3;
}

/**
 * Calculates the 3-loop beta function of g3.
 *
 * @return 3-loop beta function
 */
double SM_susy_parameters::calc_beta_g3_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g3;

   beta_g3 = Re(0.008333333333333333*threeLoop*Cube(g3)*(-523*Quad(g1) + 1635*
      Quad(g2) + 3900*Quad(g3) + 1800*Quad(Yu(2,2)) - 9*Sqr(g1)*Sqr(g2) + 616*
      Sqr(g1)*Sqr(g3) + 2520*Sqr(g2)*Sqr(g3) - 303*Sqr(g1)*Sqr(Yu(2,2)) - 1395*
      Sqr(g2)*Sqr(Yu(2,2)) - 4800*Sqr(g3)*Sqr(Yu(2,2))));


   return beta_g3;
}

/**
 * Calculates the 4-loop beta function of g3.
 *
 * @return 4-loop beta function
 */
double SM_susy_parameters::calc_beta_g3_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g3;

   beta_g3 = Re(-2472.2837425797156*Cube(g3)*Quad(oneOver16PiSqr)*(1.*Power6(g3
      ) + 0.04569149546770327*Power6(Yu(2,2)) + 0.0060672647486441755*Lambdax*
      Quad(Yu(2,2)) - 0.12603833934188147*Quad(Yu(2,2))*Sqr(g3) +
      0.16927060578749137*Quad(g3)*Sqr(Yu(2,2)) - 0.003640358849186505*Sqr(
      Lambdax)*Sqr(Yu(2,2))));


   return beta_g3;
}

/**
 * Calculates the 5-loop beta function of g3.
 *
 * @return 5-loop beta function
 */
double SM_susy_parameters::calc_beta_g3_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g3;

   beta_g3 = Re(271.4283824198132*Power11(g3)*Power5(oneOver16PiSqr));


   return beta_g3;
}

} // namespace flexiblesusy
