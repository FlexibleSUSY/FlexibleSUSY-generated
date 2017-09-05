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

// File generated at Tue 5 Sep 2017 10:39:59

#include "SM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambdax.
 *
 * @return one-loop beta function
 */
double SM_susy_parameters::calc_beta_Lambdax_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(oneOver16PiSqr*(0.27*Power(g1,4) + 2.25*Power(g2,4)
      - 9*Lambdax*Sqr(g2) + 0.9*Sqr(g1)*(-2*Lambdax + Sqr(g2)) + 4*(-3*
      traceYdAdjYdYdAdjYd - traceYeAdjYeYeAdjYe - 3*traceYuAdjYuYuAdjYu + 3*
      traceYdAdjYd*Lambdax + traceYeAdjYe*Lambdax + 3*traceYuAdjYu*Lambdax + 3*
      Sqr(Lambdax))));


   return beta_Lambdax;
}

/**
 * Calculates the two-loop beta function of Lambdax.
 *
 * @return two-loop beta function
 */
double SM_susy_parameters::calc_beta_Lambdax_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYdYdAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe =
      TRACE_STRUCT.traceYeAdjYeYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYuYuAdjYu =
      TRACE_STRUCT.traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(twoLoop*(-3.411*Power(g1,6) + 38.125*Power(g2,6) +
      60*traceYdAdjYdYdAdjYdYdAdjYd - 24*traceYdAdjYdYdAdjYuYuAdjYd + 12*
      traceYdAdjYuYuAdjYdYdAdjYd - 12*traceYdAdjYuYuAdjYuYuAdjYd + 20*
      traceYeAdjYeYeAdjYeYeAdjYe + 60*traceYuAdjYuYuAdjYuYuAdjYu - 3*
      traceYdAdjYdYdAdjYd*Lambdax - 42*traceYdAdjYuYuAdjYd*Lambdax -
      traceYeAdjYeYeAdjYe*Lambdax - 3*traceYuAdjYuYuAdjYu*Lambdax - 78*Power(
      Lambdax,3) - 0.125*Power(g2,4)*(36*traceYdAdjYd + 12*traceYeAdjYe + 36*
      traceYuAdjYu + 73*Lambdax) + 1.5*Lambdax*(15*traceYdAdjYd + 5*
      traceYeAdjYe + 15*traceYuAdjYu + 36*Lambdax)*Sqr(g2) - 0.015*Power(g1,4)*
      (-60*traceYdAdjYd + 300*traceYeAdjYe + 228*traceYuAdjYu - 629*Lambdax +
      559*Sqr(g2)) - 64*traceYdAdjYdYdAdjYd*Sqr(g3) - 64*traceYuAdjYuYuAdjYu*
      Sqr(g3) + 80*traceYdAdjYd*Lambdax*Sqr(g3) + 80*traceYuAdjYu*Lambdax*Sqr(
      g3) - 72*traceYdAdjYd*Sqr(Lambdax) - 24*traceYeAdjYe*Sqr(Lambdax) - 72*
      traceYuAdjYu*Sqr(Lambdax) - 0.025*Sqr(g1)*(289*Power(g2,4) - 6*(36*
      traceYdAdjYd + 44*traceYeAdjYe + 84*traceYuAdjYu + 39*Lambdax)*Sqr(g2) -
      4*(16*traceYdAdjYdYdAdjYd - 48*traceYeAdjYeYeAdjYe - 32*
      traceYuAdjYuYuAdjYu + 25*traceYdAdjYd*Lambdax + 75*traceYeAdjYe*Lambdax +
      85*traceYuAdjYu*Lambdax + 108*Sqr(Lambdax)))));


   return beta_Lambdax;
}

/**
 * Calculates the three-loop beta function of Lambdax.
 *
 * @return three-loop beta function
 */
double SM_susy_parameters::calc_beta_Lambdax_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = Re(0.0001*threeLoop*(-60320*Power(g1,8) - 4563640*Power
      (g2,8) - 40*Power(g1,6)*(-14084*Lambdax + 1543*Sqr(g2) - 663*Sqr(g3) -
      11117*Sqr(Yu(2,2))) + 20*Power(g2,6)*(865483*Lambdax + 15072*Sqr(g3) +
      125000*Sqr(Yu(2,2))) + 80*Power(g2,4)*(-98785*Sqr(Lambdax) - 79916*
      Lambdax*Sqr(Yu(2,2)) + Sqr(g3)*(-14286*Lambdax + 8232*Sqr(Yu(2,2))) +
      7942*Power(Yu(2,2),4)) + 2*Power(g1,4)*(130000*Power(g2,4) - 927660*Sqr(
      Lambdax) - 748599*Lambdax*Sqr(Yu(2,2)) + Sqr(g3)*(-83810*Lambdax + 20320*
      Sqr(Yu(2,2))) + 10*Sqr(g2)*(61753*Lambdax + 2210*Sqr(g3) + 21254*Sqr(Yu(2
      ,2))) + 318960*Power(Yu(2,2),4)) + 2*Sqr(g2)*(-968630*Power(Lambdax,3) +
      151443*Lambdax*Sqr(g3)*Sqr(Yu(2,2)) - 1797695*Sqr(Lambdax)*Sqr(Yu(2,2)) -
      54700*Lambdax*Power(Yu(2,2),4) + 266980*Sqr(g3)*Power(Yu(2,2),4) +
      1482760*Power(Yu(2,2),6)) - 10*Sqr(g1)*(151556*Power(g2,6) + 38745*Power(
      Lambdax,3) + 63869*Sqr(Lambdax)*Sqr(Yu(2,2)) - 4*Power(g2,4)*(39819*
      Lambdax + 1507*Sqr(g3) + 13041*Sqr(Yu(2,2))) + 42030*Lambdax*Power(Yu(2,2
      ),4) - 135720*Power(Yu(2,2),6) - 4*Sqr(g3)*(8727*Lambdax*Sqr(Yu(2,2)) +
      17570*Power(Yu(2,2),4)) + 2*Sqr(g2)*(158320*Sqr(Lambdax) - 5615*Lambdax*
      Sqr(Yu(2,2)) - 22772*Sqr(g3)*Sqr(Yu(2,2)) + 140712*Power(Yu(2,2),4))) - 5
      *(-3005675*Power(Lambdax,4) - 873000*Power(Lambdax,3)*Sqr(Yu(2,2)) -
      3536520*Sqr(Lambdax)*Power(Yu(2,2),4) + 893528*Lambdax*Power(Yu(2,2),6) +
      1945192*Power(Yu(2,2),8) + 8*Power(g3,4)*(-178484*Lambdax*Sqr(Yu(2,2)) +
      50201*Power(Yu(2,2),4)) - 4*Sqr(g3)*(80385*Sqr(Lambdax)*Sqr(Yu(2,2)) -
      662866*Lambdax*Power(Yu(2,2),4) + 500988*Power(Yu(2,2),6)))));


   return beta_Lambdax;
}

} // namespace flexiblesusy
