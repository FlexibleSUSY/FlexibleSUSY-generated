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

// File generated at Tue 22 Jan 2019 16:46:33

#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Lambdax.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_Lambdax_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_Lambdax;

   beta_Lambdax = Re(0.01*oneOver16PiSqr*(-400*g2d*g2u*gYd*gYu - 1200*
      traceYdAdjYdYdAdjYd - 400*traceYeAdjYeYeAdjYe - 1200*traceYuAdjYuYuAdjYu
      + 1200*traceYdAdjYd*Lambdax + 400*traceYeAdjYe*Lambdax + 1200*
      traceYuAdjYu*Lambdax + 27*Quad(g1) + 225*Quad(g2) - 500*Quad(g2d) - 500*
      Quad(g2u) - 100*Quad(gYd) - 100*Quad(gYu) - 180*Lambdax*Sqr(g1) - 900*
      Lambdax*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) + 600*Lambdax*Sqr(g2d) + 600*Lambdax
      *Sqr(g2u) - 200*Sqr(g2d)*Sqr(g2u) + 200*Lambdax*Sqr(gYd) - 200*Sqr(g2d)*
      Sqr(gYd) + 200*Lambdax*Sqr(gYu) - 200*Sqr(g2u)*Sqr(gYu) - 200*Sqr(gYd)*
      Sqr(gYu) + 1200*Sqr(Lambdax)));


   return beta_Lambdax;
}

/**
 * Calculates the 2-loop beta function of Lambdax.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_Lambdax_2_loop(const Susy_traces& susy_traces) const
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

   const double beta_Lambdax_1 = Re(-0.001*twoLoop*(-60000*
      traceYdAdjYdYdAdjYdYdAdjYd + 24000*traceYdAdjYdYdAdjYuYuAdjYd - 12000*
      traceYdAdjYuYuAdjYdYdAdjYd + 12000*traceYdAdjYuYuAdjYuYuAdjYd - 20000*
      traceYeAdjYeYeAdjYeYeAdjYe - 60000*traceYuAdjYuYuAdjYuYuAdjYu - 19000*g2u
      *gYd*gYu*Cube(g2d) - 19000*g2d*gYd*gYu*Cube(g2u) - 21000*g2d*g2u*gYu*Cube
      (gYd) - 21000*g2d*g2u*gYd*Cube(gYu) - 20000*g2d*g2u*gYd*gYu*Lambdax +
      3000*traceYdAdjYdYdAdjYd*Lambdax + 42000*traceYdAdjYuYuAdjYd*Lambdax +
      1000*traceYeAdjYeYeAdjYe*Lambdax + 3000*traceYuAdjYuYuAdjYu*Lambdax +
      3699*Power6(g1) - 26125*Power6(g2) - 23500*Power6(g2d) - 23500*Power6(g2u
      ) - 2500*Power6(gYd) - 2500*Power6(gYu) - 900*traceYdAdjYd*Quad(g1) +
      4500*traceYeAdjYe*Quad(g1) + 3420*traceYuAdjYu*Quad(g1) - 10035*Lambdax*
      Quad(g1) + 4500*traceYdAdjYd*Quad(g2) + 1500*traceYeAdjYe*Quad(g2) + 4500
      *traceYuAdjYu*Quad(g2) - 5875*Lambdax*Quad(g2) + 1250*Lambdax*Quad(g2d) +
      1250*Lambdax*Quad(g2u) + 250*Lambdax*Quad(gYd) + 250*Lambdax*Quad(gYu) -
      1600*traceYdAdjYdYdAdjYd*Sqr(g1) + 4800*traceYeAdjYeYeAdjYe*Sqr(g1) +
      3200*traceYuAdjYuYuAdjYu*Sqr(g1) - 2500*traceYdAdjYd*Lambdax*Sqr(g1) -
      7500*traceYeAdjYe*Lambdax*Sqr(g1) - 8500*traceYuAdjYu*Lambdax*Sqr(g1) +
      9625*Quad(g2)*Sqr(g1) + 8000*g2d*g2u*gYd*gYu*Sqr(g2) - 22500*traceYdAdjYd
      *Lambdax*Sqr(g2) - 7500*traceYeAdjYe*Lambdax*Sqr(g2) - 22500*traceYuAdjYu
      *Lambdax*Sqr(g2) + 8865*Quad(g1)*Sqr(g2) + 20000*Quad(g2d)*Sqr(g2) +
      20000*Quad(g2u)*Sqr(g2) - 5400*traceYdAdjYd*Sqr(g1)*Sqr(g2) - 6600*
      traceYeAdjYe*Sqr(g1)*Sqr(g2) - 12600*traceYuAdjYu*Sqr(g1)*Sqr(g2) - 5850*
      Lambdax*Sqr(g1)*Sqr(g2) + 270*Quad(g1)*Sqr(g2d) + 38250*Quad(g2)*Sqr(g2d)
      - 3500*Quad(g2u)*Sqr(g2d) - 8500*Quad(gYd)*Sqr(g2d) - 2250*Lambdax*Sqr(g1
      )*Sqr(g2d) - 41250*Lambdax*Sqr(g2)*Sqr(g2d) - 6300*Sqr(g1)*Sqr(g2)*Sqr(
      g2d) + 270*Quad(g1)*Sqr(g2u) + 38250*Quad(g2)*Sqr(g2u) - 3500*Quad(g2d)*
      Sqr(g2u) - 8500*Quad(gYu)*Sqr(g2u) - 2250*Lambdax*Sqr(g1)*Sqr(g2u) -
      41250*Lambdax*Sqr(g2)*Sqr(g2u) - 6300*Sqr(g1)*Sqr(g2)*Sqr(g2u) + 11000*
      Lambdax*Sqr(g2d)*Sqr(g2u) + 8000*Sqr(g2)*Sqr(g2d)*Sqr(g2u) + 64000*
      traceYdAdjYdYdAdjYd*Sqr(g3) + 64000*traceYuAdjYuYuAdjYu*Sqr(g3) - 80000*
      traceYdAdjYd*Lambdax*Sqr(g3) - 80000*traceYuAdjYu*Lambdax*Sqr(g3) + 90*
      Quad(g1)*Sqr(gYd) + 750*Quad(g2)*Sqr(gYd) - 5500*Quad(g2d)*Sqr(gYd) -
      8500*Quad(gYu)*Sqr(gYd) - 750*Lambdax*Sqr(g1)*Sqr(gYd) - 3750*Lambdax*Sqr
      (g2)*Sqr(gYd) + 300*Sqr(g1)*Sqr(g2)*Sqr(gYd) + 500*Lambdax*Sqr(g2d)*Sqr(
      gYd) + 4000*Sqr(g2)*Sqr(g2d)*Sqr(gYd) - 10500*Sqr(g2d)*Sqr(g2u)*Sqr(gYd)
      + 90*Quad(g1)*Sqr(gYu) + 750*Quad(g2)*Sqr(gYu) - 5500*Quad(g2u)*Sqr(gYu)
      - 8500*Quad(gYd)*Sqr(gYu) - 750*Lambdax*Sqr(g1)*Sqr(gYu) - 3750*Lambdax*
      Sqr(g2)*Sqr(gYu) + 300*Sqr(g1)*Sqr(g2)*Sqr(gYu) + 500*Lambdax*Sqr(g2u)*
      Sqr(gYu) + 4000*Sqr(g2)*Sqr(g2u)*Sqr(gYu) - 10500*Sqr(g2d)*Sqr(g2u)*Sqr(
      gYu) - 3000*Lambdax*Sqr(gYd)*Sqr(gYu) - 9500*Sqr(g2d)*Sqr(gYd)*Sqr(gYu) -
      9500*Sqr(g2u)*Sqr(gYd)*Sqr(gYu) - 10800*Sqr(g1)*Sqr(Lambdax) - 54000*Sqr(
      g2)*Sqr(Lambdax) + 36000*Sqr(g2d)*Sqr(Lambdax) + 36000*Sqr(g2u)*Sqr(
      Lambdax) + 12000*Sqr(gYd)*Sqr(Lambdax) + 12000*Sqr(gYu)*Sqr(Lambdax)));
   const double beta_Lambdax_2 = Re(-6*twoLoop*(12*traceYdAdjYd + 4*
      traceYeAdjYe + 12*traceYuAdjYu + 13*Lambdax)*Sqr(Lambdax));

   beta_Lambdax = beta_Lambdax_1 + beta_Lambdax_2;


   return beta_Lambdax;
}

/**
 * Calculates the 3-loop beta function of Lambdax.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_Lambdax_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

/**
 * Calculates the 4-loop beta function of Lambdax.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_Lambdax_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

/**
 * Calculates the 5-loop beta function of Lambdax.
 *
 * @return 5-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_Lambdax_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Lambdax;

   beta_Lambdax = 0;


   return beta_Lambdax;
}

} // namespace flexiblesusy
