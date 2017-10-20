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

// File generated at Fri 20 Oct 2017 08:35:10

#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambda1.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_Lambda1;

   beta_Lambda1 = Re(oneOver16PiSqr*(2*Lambda3*Lambda4 + 12*Lambda1*
      traceYdAdjYd - 6*traceYdAdjYdYdAdjYd + 4*Lambda1*traceYeAdjYe - 2*
      traceYeAdjYeYeAdjYe + 0.135*Quad(g1) - 2.5*Quad(g1d) - 0.5*Quad(g1dp) +
      1.125*Quad(g2) + 2*Lambda1*Sqr(g1dp) - Sqr(g1d)*(-6*Lambda1 + Sqr(g1dp))
      - 9*Lambda1*Sqr(g2) + 0.45*Sqr(g1)*(-4*Lambda1 + Sqr(g2)) + 24*Sqr(
      Lambda1) + 2*Sqr(Lambda3) + Sqr(Lambda4) + Sqr(Lambda5) + 12*Sqr(Lambda6)
      ));


   return beta_Lambda1;
}

/**
 * Calculates the 2-loop beta function of Lambda1.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYdAdjYdYdAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd =
      TRACE_STRUCT.traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe =
      TRACE_STRUCT.traceYeAdjYeYeAdjYeYeAdjYe;


   double beta_Lambda1;

   const double beta_Lambda1_1 = Re(-0.0025*twoLoop*(765*Power6(g1) + 3*
      Quad(g1)*(-1382*Lambda1 - 120*Lambda3 - 60*Lambda4 - 60*traceYdAdjYd +
      300*traceYeAdjYe + 18*Sqr(g1d) + 6*Sqr(g1dp) + 605*Sqr(g2)) - 5*Sqr(g1)*(
      192*Lambda3*Lambda4 + 200*Lambda1*traceYdAdjYd + 64*traceYdAdjYdYdAdjYd -
      399*Quad(g2) + 468*Lambda1*Sqr(g2) + 120*Lambda4*Sqr(g2) + 216*
      traceYdAdjYd*Sqr(g2) + 264*traceYeAdjYe*Sqr(g2) - 12*Sqr(g1dp)*(-5*
      Lambda1 + Sqr(g2)) + 36*Sqr(g1d)*(5*Lambda1 + 7*Sqr(g2)) + 1728*Sqr(
      Lambda1) + 192*Sqr(Lambda3) + 96*Sqr(Lambda4) - 48*Sqr(Lambda5) + 864*Sqr
      (Lambda6)) - 25*(-320*Lambda1*Lambda3*Lambda4 - 576*Lambda3*Lambda6*
      Lambda7 - 448*Lambda4*Lambda6*Lambda7 - 320*Lambda5*Lambda6*Lambda7 - 48*
      Lambda1*traceYdAdjYdYdAdjYd + 480*traceYdAdjYdYdAdjYdYdAdjYd - 144*
      Lambda1*traceYdAdjYuYuAdjYd + 96*traceYdAdjYuYuAdjYdYdAdjYd - 4992*Cube(
      Lambda1) - 128*Cube(Lambda3) - 96*Cube(Lambda4) + 188*Power6(g1d) + 20*
      Power6(g1dp) + 195*Power6(g2) + 138*Lambda1*Quad(g2) + 120*Lambda3*Quad(
      g2) + 60*Lambda4*Quad(g2) - 36*traceYdAdjYd*Quad(g2) + 192*Lambda3*
      Lambda4*Sqr(g2) + 360*Lambda1*traceYdAdjYd*Sqr(g2) + 4*Quad(g1d)*(11*Sqr(
      g1dp) - 5*(Lambda1 + 8*Sqr(g2) - 2*Sqr(g2u))) - 96*Lambda3*Lambda4*Sqr(
      g2u) - 32*Lambda3*Lambda4*Sqr(g2up) + Quad(g1dp)*(-4*Lambda1 + 8*Sqr(g2up
      )) + 1280*Lambda1*traceYdAdjYd*Sqr(g3) - 512*traceYdAdjYdYdAdjYd*Sqr(g3)
      - 2304*traceYdAdjYd*Sqr(Lambda1) + 1728*Sqr(g2)*Sqr(Lambda1) - 320*
      Lambda1*Sqr(Lambda3) - 192*Lambda4*Sqr(Lambda3) + 192*Sqr(g2)*Sqr(Lambda3
      ) - 96*Sqr(g2u)*Sqr(Lambda3) - 32*Sqr(g2up)*Sqr(Lambda3) - 192*Lambda1*
      Sqr(Lambda4) - 256*Lambda3*Sqr(Lambda4) + 48*Sqr(g2)*Sqr(Lambda4) - 48*
      Sqr(g2u)*Sqr(Lambda4) - 16*Sqr(g2up)*Sqr(Lambda4) - 224*Lambda1*Sqr(
      Lambda5) - 320*Lambda3*Sqr(Lambda5) - 352*Lambda4*Sqr(Lambda5) - 48*Sqr(
      g2u)*Sqr(Lambda5) - 16*Sqr(g2up)*Sqr(Lambda5) - 5088*Lambda1*Sqr(Lambda6)
      - 1056*Lambda3*Sqr(Lambda6) - 1120*Lambda4*Sqr(Lambda6) - 1184*Lambda5*
      Sqr(Lambda6) - 576*traceYdAdjYd*Sqr(Lambda6) + 864*Sqr(g2)*Sqr(Lambda6) -
      288*Sqr(g2u)*Sqr(Lambda6) - 96*Sqr(g2up)*Sqr(Lambda6) - 6*Sqr(g1dp)*(
      Quad(g2) - 10*Lambda1*Sqr(g2) + 4*Lambda1*Sqr(g2up) + 64*Sqr(Lambda1) +
      16*Sqr(Lambda6)) + 2*Sqr(g1d)*(34*Quad(g1dp) - 4*Sqr(g1dp)*(Lambda1 + 4*
      Sqr(g2) - Sqr(g2u) - Sqr(g2up)) - 3*(51*Quad(g2) - 110*Lambda1*Sqr(g2) +
      12*Lambda1*Sqr(g2u) + 192*Sqr(Lambda1) + 48*Sqr(Lambda6))) + 96*Lambda1*
      Sqr(Lambda7) - 288*Lambda3*Sqr(Lambda7) - 224*Lambda4*Sqr(Lambda7) - 160*
      Lambda5*Sqr(Lambda7))));
   const double beta_Lambda1_2 = Re(-0.05*twoLoop*(15*traceYeAdjYe*Quad(
      g2) - 150*Lambda1*traceYeAdjYe*Sqr(g2) + 2*((-75*Lambda1*traceYeAdjYe +
      24*traceYeAdjYeYeAdjYe)*Sqr(g1) + 10*(Lambda1*traceYeAdjYeYeAdjYe + 48*
      traceYeAdjYe*Sqr(Lambda1) + 2*(-5*traceYeAdjYeYeAdjYeYeAdjYe + 3*
      traceYuAdjYu*(2*Lambda3*Lambda4 + 2*Sqr(Lambda3) + Sqr(Lambda4) + Sqr(
      Lambda5)) + 6*(traceYeAdjYe + 3*traceYuAdjYu)*Sqr(Lambda6))))));

   beta_Lambda1 = beta_Lambda1_1 + beta_Lambda1_2;


   return beta_Lambda1;
}

/**
 * Calculates the 3-loop beta function of Lambda1.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_Lambda1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambda1;

   beta_Lambda1 = 0;


   return beta_Lambda1;
}

} // namespace flexiblesusy
