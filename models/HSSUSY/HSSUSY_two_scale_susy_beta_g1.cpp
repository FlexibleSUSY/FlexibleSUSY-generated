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

// File generated at Tue 5 Sep 2017 10:38:53

#include "HSSUSY_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g1.
 *
 * @return one-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g1_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g1;

   beta_g1 = Re(4.1*Power(g1,3)*oneOver16PiSqr);


   return beta_g1;
}

/**
 * Calculates the two-loop beta function of g1.
 *
 * @return two-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g1_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g1;

   beta_g1 = Re(0.02*Power(g1,3)*twoLoop*(199*Sqr(g1) + 5*(-5*
      traceYdAdjYd - 15*traceYeAdjYe - 17*traceYuAdjYu + 27*Sqr(g2) + 88*Sqr(g3
      ))));


   return beta_g1;
}

/**
 * Calculates the three-loop beta function of g1.
 *
 * @return three-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g1_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1;

   beta_g1 = Re(-0.000041666666666666665*Power(g1,3)*threeLoop*(388613*
      Power(g1,4) - 10*Sqr(g1)*(648*Lambdax + 1845*Sqr(g2) - 4384*Sqr(g3) -
      8481*Sqr(Yu(2,2))) - 75*(3945*Power(g2,4) - 6*Sqr(g2)*(-24*Lambdax + 32*
      Sqr(g3) + 785*Sqr(Yu(2,2))) + 4*(4752*Power(g3,4) - 36*Sqr(Lambdax) - 464
      *Sqr(g3)*Sqr(Yu(2,2)) + 945*Power(Yu(2,2),4)))));


   return beta_g1;
}

} // namespace flexiblesusy
