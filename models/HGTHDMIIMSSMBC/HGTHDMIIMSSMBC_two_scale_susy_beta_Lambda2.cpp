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

// File generated at Mon 9 May 2016 12:04:16

#include "HGTHDMIIMSSMBC_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Lambda2.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_one_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambda2;

   beta_Lambda2 = Re(oneOver16PiSqr*(0.135*Power(g1,4) + 1.125*Power(g2,4
      ) - 2.5*Power(g2u,4) - 0.5*Power(g2up,4) + 2*Lambda3*Lambda4 + 12*Lambda2
      *traceYuAdjYu - 6*traceYuAdjYuYuAdjYu - 1.8*Lambda2*Sqr(g1) - 9*Lambda2*
      Sqr(g2) + 0.45*Sqr(g1)*Sqr(g2) + 6*Lambda2*Sqr(g2u) + 2*Lambda2*Sqr(g2up)
      - Sqr(g2u)*Sqr(g2up) + 24*Sqr(Lambda2) + 2*Sqr(Lambda3) + Sqr(Lambda4) +
      Sqr(Lambda5) + 12*Sqr(Lambda7)));


   return beta_Lambda2;
}

/**
 * Calculates the two-loop beta function of Lambda2.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_two_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYuYuAdjYuYuAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYuYuAdjYu =
      TRACE_STRUCT.traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_Lambda2;

   const double beta_Lambda2_1 = Re(-0.0025*twoLoop*(765*Power(g1,6) + 3*
      Power(g1,4)*(605*Sqr(g2) + 2*(-691*Lambda2 - 60*Lambda3 - 30*Lambda4 +
      114*traceYuAdjYu + 9*Sqr(g2u) + 3*Sqr(g2up))) + 5*Sqr(g1)*(399*Power(g2,4
      ) - 12*Sqr(g2)*(39*Lambda2 + 10*Lambda4 + 42*traceYuAdjYu + 21*Sqr(g2u) -
      Sqr(g2up)) - 4*(48*Lambda3*Lambda4 + 170*Lambda2*traceYuAdjYu + 45*
      Lambda2*Sqr(g2u) + 15*Lambda2*Sqr(g2up) + 432*Sqr(Lambda2) + 48*Sqr(
      Lambda3) + 24*Sqr(Lambda4) - 12*Sqr(Lambda5) + 216*Sqr(Lambda7))) - 25*(
      195*Power(g2,6) - 6*Power(g2,4)*(-23*Lambda2 - 20*Lambda3 - 10*Lambda4 +
      6*traceYuAdjYu + 51*Sqr(g2u) + Sqr(g2up)) + 4*(47*Power(g2u,6) + 5*Power(
      g2up,6) - Power(g2up,4)*Lambda2 - 1248*Power(Lambda2,3) - 32*Power(
      Lambda3,3) - 80*Lambda2*Lambda3*Lambda4 - 24*Power(Lambda4,3) - 144*
      Lambda3*Lambda6*Lambda7 - 112*Lambda4*Lambda6*Lambda7 - 80*Lambda5*
      Lambda6*Lambda7 - 48*Lambda3*Lambda4*traceYdAdjYd - 36*Lambda2*
      traceYdAdjYuYuAdjYd + 24*traceYdAdjYuYuAdjYuYuAdjYd - 16*Lambda3*Lambda4*
      traceYeAdjYe + 2*Power(g2up,4)*Sqr(g1dp) - 8*Lambda3*Lambda4*Sqr(g1dp) -
      6*Lambda2*Sqr(g1dp)*Sqr(g2up) + Power(g2u,4)*(-5*Lambda2 + 11*Sqr(g2up))
      - 96*Sqr(g2up)*Sqr(Lambda2) - 80*Lambda2*Sqr(Lambda3) - 48*Lambda4*Sqr(
      Lambda3) - 48*traceYdAdjYd*Sqr(Lambda3) - 16*traceYeAdjYe*Sqr(Lambda3) -
      8*Sqr(g1dp)*Sqr(Lambda3) - 48*Lambda2*Sqr(Lambda4) - 64*Lambda3*Sqr(
      Lambda4) - 24*traceYdAdjYd*Sqr(Lambda4) - 8*traceYeAdjYe*Sqr(Lambda4) - 4
      *Sqr(g1dp)*Sqr(Lambda4) - 56*Lambda2*Sqr(Lambda5) - 80*Lambda3*Sqr(
      Lambda5) - 88*Lambda4*Sqr(Lambda5) - 24*traceYdAdjYd*Sqr(Lambda5) - 8*
      traceYeAdjYe*Sqr(Lambda5) - 4*Sqr(g1dp)*Sqr(Lambda5) + 24*Lambda2*Sqr(
      Lambda6) - 72*Lambda3*Sqr(Lambda6) - 56*Lambda4*Sqr(Lambda6) - 40*Lambda5
      *Sqr(Lambda6) - 1272*Lambda2*Sqr(Lambda7) - 264*Lambda3*Sqr(Lambda7) -
      280*Lambda4*Sqr(Lambda7) - 296*Lambda5*Sqr(Lambda7) - 144*traceYdAdjYd*
      Sqr(Lambda7) - 48*traceYeAdjYe*Sqr(Lambda7) - 24*Sqr(g1dp)*Sqr(Lambda7) -
      24*Sqr(g2up)*Sqr(Lambda7) + Sqr(g2u)*(17*Power(g2up,4) - 2*Lambda2*Sqr(
      g2up) + 2*Sqr(g1dp)*Sqr(g2up) - 72*(4*Sqr(Lambda2) + Sqr(Lambda7))) + 2*
      Sqr(g1d)*(5*Power(g2u,4) + Sqr(g2u)*(-9*Lambda2 + Sqr(g2up)) - 6*(2*
      Lambda3*Lambda4 + 2*Sqr(Lambda3) + Sqr(Lambda4) + Sqr(Lambda5) + 6*Sqr(
      Lambda7)))) + 4*Sqr(g2)*(-40*Power(g2u,4) + Sqr(g2u)*(165*Lambda2 - 8*Sqr
      (g2up)) + 3*(5*Lambda2*Sqr(g2up) + 4*(4*Lambda3*Lambda4 + 36*Sqr(Lambda2)
      + 4*Sqr(Lambda3) + Sqr(Lambda4) + 18*Sqr(Lambda7)))))));
   const double beta_Lambda2_2 = Re(twoLoop*(-3*Lambda2*
      traceYuAdjYuYuAdjYu + 30*traceYuAdjYuYuAdjYuYuAdjYu - 1.6*
      traceYuAdjYuYuAdjYu*Sqr(g1) + 22.5*Lambda2*traceYuAdjYu*Sqr(g2) + 16*(5*
      Lambda2*traceYuAdjYu - 2*traceYuAdjYuYuAdjYu)*Sqr(g3) - 144*traceYuAdjYu*
      Sqr(Lambda2) - 36*traceYuAdjYu*Sqr(Lambda7)));

   beta_Lambda2 = beta_Lambda2_1 + beta_Lambda2_2;


   return beta_Lambda2;
}

/**
 * Calculates the three-loop beta function of Lambda2.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda2_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda2;

   beta_Lambda2 = 0;


   return beta_Lambda2;
}

} // namespace flexiblesusy
