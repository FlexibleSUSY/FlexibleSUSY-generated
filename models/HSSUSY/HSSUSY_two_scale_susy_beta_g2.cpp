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

// File generated at Tue 5 Sep 2017 10:38:54

#include "HSSUSY_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g2.
 *
 * @return one-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g2_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = Re(-3.1666666666666665*Power(g2,3)*oneOver16PiSqr);


   return beta_g2;
}

/**
 * Calculates the two-loop beta function of g2.
 *
 * @return two-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g2_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2;

   beta_g2 = Re(0.03333333333333333*Power(g2,3)*twoLoop*(27*Sqr(g1) + 5*(
      -3*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu) + 35*Sqr(g2) + 72*Sqr
      (g3))));


   return beta_g2;
}

/**
 * Calculates the three-loop beta function of g2.
 *
 * @return three-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_g2_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = Re(0.000023148148148148147*Power(g2,3)*threeLoop*(-151119*
      Power(g1,4) + 270*Sqr(g1)*(24*Lambdax + 873*Sqr(g2) - 32*Sqr(g3) - 593*
      Sqr(Yu(2,2))) + 25*(324953*Power(g2,4) + 162*Sqr(g2)*(8*Lambdax + 416*Sqr
      (g3) - 243*Sqr(Yu(2,2))) + 108*(1296*Power(g3,4) - 12*Sqr(Lambdax) - 112*
      Sqr(g3)*Sqr(Yu(2,2)) + 147*Power(Yu(2,2),4)))));


   return beta_g2;
}

} // namespace flexiblesusy
