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


#include "SplitMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of gYu.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_gYu;

   beta_gYu = Re(0.05*(60*g2d*g2u*gYd + 60*gYu*traceYdAdjYd + 20*gYu*
      traceYeAdjYe + 60*gYu*traceYuAdjYu + 25*Cube(gYu) - 9*gYu*Sqr(g1) - 45*
      gYu*Sqr(g2) + 30*gYu*Sqr(g2d) + 45*gYu*Sqr(g2u) + 40*gYu*Sqr(gYd)));


   return oneLoop * beta_gYu;
}

/**
 * Calculates the 2-loop beta function of gYu.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_gYu;

   beta_gYu = Re(0.00125*(-7200*g2d*g2u*gYd*traceYdAdjYd - 5400*gYu*
      traceYdAdjYdYdAdjYd + 1200*gYu*traceYdAdjYuYuAdjYd - 2400*g2d*g2u*gYd*
      traceYeAdjYe - 1800*gYu*traceYeAdjYeYeAdjYe - 7200*g2d*g2u*gYd*
      traceYuAdjYu - 5400*gYu*traceYuAdjYuYuAdjYu - 6600*g2u*gYd*Cube(g2d) -
      3600*g2d*gYd*Cube(g2u) - 1800*g2d*g2u*Cube(gYd) - 2700*traceYdAdjYd*Cube(
      gYu) - 900*traceYeAdjYe*Cube(gYu) - 2700*traceYuAdjYu*Cube(gYu) - 2400*
      g2d*g2u*gYd*Lambdax - 2400*Cube(gYu)*Lambdax - 600*Power5(gYu) + 468*gYu*
      Quad(g1) - 3400*gYu*Quad(g2) - 2250*gYu*Quad(g2d) - 4950*gYu*Quad(g2u) -
      1800*gYu*Quad(gYd) + 360*g2d*g2u*gYd*Sqr(g1) + 500*gYu*traceYdAdjYd*Sqr(
      g1) + 1500*gYu*traceYeAdjYe*Sqr(g1) + 1700*gYu*traceYuAdjYu*Sqr(g1) +
      1545*Cube(gYu)*Sqr(g1) + 10200*g2d*g2u*gYd*Sqr(g2) + 4500*gYu*
      traceYdAdjYd*Sqr(g2) + 1500*gYu*traceYeAdjYe*Sqr(g2) + 4500*gYu*
      traceYuAdjYu*Sqr(g2) + 4125*Cube(gYu)*Sqr(g2) - 1080*gYu*Sqr(g1)*Sqr(g2)
      - 1350*Cube(gYu)*Sqr(g2d) + 450*gYu*Sqr(g1)*Sqr(g2d) + 8250*gYu*Sqr(g2)*
      Sqr(g2d) - 2700*gYu*traceYdAdjYd*Sqr(g2u) - 900*gYu*traceYeAdjYe*Sqr(g2u)
      - 2700*gYu*traceYuAdjYu*Sqr(g2u) - 450*Cube(gYu)*Sqr(g2u) - 2400*gYu*
      Lambdax*Sqr(g2u) + 945*gYu*Sqr(g1)*Sqr(g2u) + 13725*gYu*Sqr(g2)*Sqr(g2u)
      - 2100*gYu*Sqr(g2d)*Sqr(g2u) + 16000*gYu*traceYdAdjYd*Sqr(g3) + 16000*gYu
      *traceYuAdjYu*Sqr(g3) - 4200*gYu*traceYdAdjYd*Sqr(gYd) - 1400*gYu*
      traceYeAdjYe*Sqr(gYd) - 4200*gYu*traceYuAdjYu*Sqr(gYd) - 3000*Cube(gYu)*
      Sqr(gYd) - 2400*gYu*Lambdax*Sqr(gYd) + 60*gYu*Sqr(g1)*Sqr(gYd) + 3900*gYu
      *Sqr(g2)*Sqr(gYd) - 3750*gYu*Sqr(g2d)*Sqr(gYd) - 3750*gYu*Sqr(g2u)*Sqr(
      gYd) - 4800*g2d*g2u*gYd*Sqr(gYu) + 1200*gYu*Sqr(Lambdax)));


   return twoLoop * beta_gYu;
}

/**
 * Calculates the 3-loop beta function of gYu.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYu;

   beta_gYu = 0;


   return threeLoop * beta_gYu;
}

/**
 * Calculates the 4-loop beta function of gYu.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYu;

   beta_gYu = 0;


   return fourLoop * beta_gYu;
}

/**
 * Calculates the 5-loop beta function of gYu.
 *
 * @return 5-loop beta function
 */
double SplitMSSM_susy_parameters::calc_beta_gYu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gYu;

   beta_gYu = 0;


   return fiveLoop * beta_gYu;
}

} // namespace flexiblesusy
