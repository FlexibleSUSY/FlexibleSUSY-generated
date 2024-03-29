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


#include "HSSUSY_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2.
 *
 * @return 1-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g2_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = Re(-3.1666666666666665*Cube(g2));


   return oneLoop * beta_g2;
}

/**
 * Calculates the 2-loop beta function of g2.
 *
 * @return 2-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g2_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2;

   beta_g2 = Re(0.03333333333333333*Cube(g2)*(-45*traceYdAdjYd - 15*
      traceYeAdjYe - 45*traceYuAdjYu + 27*Sqr(g1) + 175*Sqr(g2) + 360*Sqr(g3)))
      ;


   return twoLoop * beta_g2;
}

/**
 * Calculates the 3-loop beta function of g2.
 *
 * @return 3-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g2_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = Re(0.000023148148148148147*Cube(g2)*(-151119*Quad(g1) + 8123825*
      Quad(g2) + 3499200*Quad(g3) + 396900*Quad(Yu(2,2)) + 6480*Lambdax*Sqr(g1)
      + 32400*Lambdax*Sqr(g2) + 235710*Sqr(g1)*Sqr(g2) - 8640*Sqr(g1)*Sqr(g3) +
      1684800*Sqr(g2)*Sqr(g3) - 32400*Sqr(Lambdax) - 160110*Sqr(g1)*Sqr(Yu(2,2)
      ) - 984150*Sqr(g2)*Sqr(Yu(2,2)) - 302400*Sqr(g3)*Sqr(Yu(2,2))));


   return threeLoop * beta_g2;
}

/**
 * Calculates the 4-loop beta function of g2.
 *
 * @return 4-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g2_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return fourLoop * beta_g2;
}

/**
 * Calculates the 5-loop beta function of g2.
 *
 * @return 5-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g2_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return fiveLoop * beta_g2;
}

} // namespace flexiblesusy
