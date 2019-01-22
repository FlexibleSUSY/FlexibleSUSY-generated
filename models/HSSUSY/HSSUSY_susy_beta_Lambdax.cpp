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

// File generated at Tue 22 Jan 2019 16:35:36

#include "HSSUSY_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambdax.
 *
 * @return 1-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(0.01*oneOver16PiSqr*(-1200*traceYdAdjYdYdAdjYd - 400*
      traceYeAdjYeYeAdjYe - 1200*traceYuAdjYuYuAdjYu + 1200*traceYdAdjYd*
      Lambdax + 400*traceYeAdjYe*Lambdax + 1200*traceYuAdjYu*Lambdax + 27*Quad(
      g1) + 225*Quad(g2) - 180*Lambdax*Sqr(g1) - 900*Lambdax*Sqr(g2) + 90*Sqr(
      g1)*Sqr(g2) + 1200*Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the 2-loop beta function of Lambdax.
 *
 * @return 2-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYdYdAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd = TRACE_STRUCT.
      traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe = TRACE_STRUCT.
      traceYeAdjYeYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYuYuAdjYu = TRACE_STRUCT.
      traceYuAdjYuYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(0.001*twoLoop*(60000*traceYdAdjYdYdAdjYdYdAdjYd - 24000*
      traceYdAdjYdYdAdjYuYuAdjYd + 12000*traceYdAdjYuYuAdjYdYdAdjYd - 12000*
      traceYdAdjYuYuAdjYuYuAdjYd + 20000*traceYeAdjYeYeAdjYeYeAdjYe + 60000*
      traceYuAdjYuYuAdjYuYuAdjYu - 78000*Cube(Lambdax) - 3000*
      traceYdAdjYdYdAdjYd*Lambdax - 42000*traceYdAdjYuYuAdjYd*Lambdax - 1000*
      traceYeAdjYeYeAdjYe*Lambdax - 3000*traceYuAdjYuYuAdjYu*Lambdax - 3411*
      Power6(g1) + 38125*Power6(g2) + 900*traceYdAdjYd*Quad(g1) - 4500*
      traceYeAdjYe*Quad(g1) - 3420*traceYuAdjYu*Quad(g1) + 9435*Lambdax*Quad(g1
      ) - 4500*traceYdAdjYd*Quad(g2) - 1500*traceYeAdjYe*Quad(g2) - 4500*
      traceYuAdjYu*Quad(g2) - 9125*Lambdax*Quad(g2) + 1600*traceYdAdjYdYdAdjYd*
      Sqr(g1) - 4800*traceYeAdjYeYeAdjYe*Sqr(g1) - 3200*traceYuAdjYuYuAdjYu*Sqr
      (g1) + 2500*traceYdAdjYd*Lambdax*Sqr(g1) + 7500*traceYeAdjYe*Lambdax*Sqr(
      g1) + 8500*traceYuAdjYu*Lambdax*Sqr(g1) - 7225*Quad(g2)*Sqr(g1) + 22500*
      traceYdAdjYd*Lambdax*Sqr(g2) + 7500*traceYeAdjYe*Lambdax*Sqr(g2) + 22500*
      traceYuAdjYu*Lambdax*Sqr(g2) - 8385*Quad(g1)*Sqr(g2) + 5400*traceYdAdjYd*
      Sqr(g1)*Sqr(g2) + 6600*traceYeAdjYe*Sqr(g1)*Sqr(g2) + 12600*traceYuAdjYu*
      Sqr(g1)*Sqr(g2) + 5850*Lambdax*Sqr(g1)*Sqr(g2) - 64000*
      traceYdAdjYdYdAdjYd*Sqr(g3) - 64000*traceYuAdjYuYuAdjYu*Sqr(g3) + 80000*
      traceYdAdjYd*Lambdax*Sqr(g3) + 80000*traceYuAdjYu*Lambdax*Sqr(g3) - 72000
      *traceYdAdjYd*Sqr(Lambdax) - 24000*traceYeAdjYe*Sqr(Lambdax) - 72000*
      traceYuAdjYu*Sqr(Lambdax) + 10800*Sqr(g1)*Sqr(Lambdax) + 54000*Sqr(g2)*
      Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the 3-loop beta function of Lambdax.
 *
 * @return 3-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = Re(0.0001*threeLoop*(563360*Lambdax*Power6(g1) + 17309660*
      Lambdax*Power6(g2) - 4467640*Lambdax*Power6(Yu(2,2)) - 60320*Power8(g1) -
      4563640*Power8(g2) - 9725960*Power8(Yu(2,2)) + 260000*Quad(g1)*Quad(g2) +
      15028375*Quad(Lambdax) + 637920*Quad(g1)*Quad(Yu(2,2)) + 635360*Quad(g2)*
      Quad(Yu(2,2)) - 2008040*Quad(g3)*Quad(Yu(2,2)) - 387450*Cube(Lambdax)*Sqr
      (g1) - 1515560*Power6(g2)*Sqr(g1) + 1357200*Power6(Yu(2,2))*Sqr(g1) +
      1592760*Lambdax*Quad(g2)*Sqr(g1) - 420300*Lambdax*Quad(Yu(2,2))*Sqr(g1) -
      1937260*Cube(Lambdax)*Sqr(g2) - 61720*Power6(g1)*Sqr(g2) + 2965520*Power6
      (Yu(2,2))*Sqr(g2) + 1235060*Lambdax*Quad(g1)*Sqr(g2) - 109400*Lambdax*
      Quad(Yu(2,2))*Sqr(g2) - 2814240*Quad(Yu(2,2))*Sqr(g1)*Sqr(g2) + 26520*
      Power6(g1)*Sqr(g3) + 301440*Power6(g2)*Sqr(g3) + 10019760*Power6(Yu(2,2))
      *Sqr(g3) - 167620*Lambdax*Quad(g1)*Sqr(g3) - 1142880*Lambdax*Quad(g2)*Sqr
      (g3) - 13257320*Lambdax*Quad(Yu(2,2))*Sqr(g3) + 60280*Quad(g2)*Sqr(g1)*
      Sqr(g3) + 702800*Quad(Yu(2,2))*Sqr(g1)*Sqr(g3) + 44200*Quad(g1)*Sqr(g2)*
      Sqr(g3) + 533960*Quad(Yu(2,2))*Sqr(g2)*Sqr(g3) - 1855320*Quad(g1)*Sqr(
      Lambdax) - 7902800*Quad(g2)*Sqr(Lambdax) + 17682600*Quad(Yu(2,2))*Sqr(
      Lambdax) - 3166400*Sqr(g1)*Sqr(g2)*Sqr(Lambdax) + 4365000*Cube(Lambdax)*
      Sqr(Yu(2,2)) + 444680*Power6(g1)*Sqr(Yu(2,2)) + 2500000*Power6(g2)*Sqr(Yu
      (2,2)) - 1497198*Lambdax*Quad(g1)*Sqr(Yu(2,2)) - 6393280*Lambdax*Quad(g2)
      *Sqr(Yu(2,2)) + 7139360*Lambdax*Quad(g3)*Sqr(Yu(2,2)) + 521640*Quad(g2)*
      Sqr(g1)*Sqr(Yu(2,2)) + 425080*Quad(g1)*Sqr(g2)*Sqr(Yu(2,2)) + 112300*
      Lambdax*Sqr(g1)*Sqr(g2)*Sqr(Yu(2,2)) + 40640*Quad(g1)*Sqr(g3)*Sqr(Yu(2,2)
      ) + 658560*Quad(g2)*Sqr(g3)*Sqr(Yu(2,2)) + 349080*Lambdax*Sqr(g1)*Sqr(g3)
      *Sqr(Yu(2,2)) + 302886*Lambdax*Sqr(g2)*Sqr(g3)*Sqr(Yu(2,2)) + 455440*Sqr(
      g1)*Sqr(g2)*Sqr(g3)*Sqr(Yu(2,2)) - 638690*Sqr(g1)*Sqr(Lambdax)*Sqr(Yu(2,2
      )) - 3595390*Sqr(g2)*Sqr(Lambdax)*Sqr(Yu(2,2)) + 1607700*Sqr(g3)*Sqr(
      Lambdax)*Sqr(Yu(2,2))));


   return beta_Lambdax;
}

/**
 * Calculates the 4-loop beta function of Lambdax.
 *
 * @return 4-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = Re(16616.34*Power6(g3)*Quad(oneOver16PiSqr)*Quad(Yu(2,2)));


   return beta_Lambdax;
}

/**
 * Calculates the 5-loop beta function of Lambdax.
 *
 * @return 5-loop beta function
 */
double HSSUSY_susy_parameters::calc_beta_Lambdax_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
