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


#include "HGTHDMIIMSSMBC_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g2up.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2up_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_g2up;

   beta_g2up = Re(0.05*g2up*(60*traceYuAdjYu - 9*Sqr(g1) + 10*Sqr(g1dp) - 45*
      Sqr(g2) + 45*Sqr(g2u) + 25*Sqr(g2up)));


   return oneLoop * beta_g2up;
}

/**
 * Calculates the 2-loop beta function of g2up.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2up_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_g2up;

   beta_g2up = Re(0.00125*(-2400*g1d*g1dp*g2u*Lambda4 + 800*g2up*Lambda3*
      Lambda4 - 1800*g2up*traceYdAdjYuYuAdjYd - 5400*g2up*traceYuAdjYuYuAdjYu +
      1200*g2up*AbsSqr(Lambda5) + 1200*g2up*AbsSqr(Lambda6) + 3600*g2up*AbsSqr(
      Lambda7) - 4800*Lambda2*Cube(g2up) - 2700*traceYuAdjYu*Cube(g2up) - 600*
      Power5(g2up) + 516*g2up*Quad(g1) - 350*g2up*Quad(g1dp) - 3000*g2up*Quad(
      g2) - 4950*g2up*Quad(g2u) + 1700*g2up*traceYuAdjYu*Sqr(g1) + 1545*Cube(
      g2up)*Sqr(g1) - 1600*g2up*Lambda3*Sqr(g1dp) - 800*g2up*Lambda4*Sqr(g1dp)
      - 1800*g2up*traceYdAdjYd*Sqr(g1dp) - 600*g2up*traceYeAdjYe*Sqr(g1dp) -
      350*Cube(g2up)*Sqr(g1dp) - 210*g2up*Sqr(g1)*Sqr(g1dp) - 1050*g2up*Sqr(g1d
      )*Sqr(g1dp) - 3600*g1d*g1dp*g2u*Sqr(g2) + 4500*g2up*traceYuAdjYu*Sqr(g2)
      + 4125*Cube(g2up)*Sqr(g2) - 1080*g2up*Sqr(g1)*Sqr(g2) + 2550*g2up*Sqr(
      g1dp)*Sqr(g2) - 4800*g2up*Lambda2*Sqr(g2u) - 2700*g2up*traceYuAdjYu*Sqr(
      g2u) - 450*Cube(g2up)*Sqr(g2u) + 945*g2up*Sqr(g1)*Sqr(g2u) - 1050*g2up*
      Sqr(g1d)*Sqr(g2u) + 13725*g2up*Sqr(g2)*Sqr(g2u) + 16000*g2up*traceYuAdjYu
      *Sqr(g3) + 4800*g2up*Sqr(Lambda2) + 800*g2up*Sqr(Lambda3) + 800*g2up*Sqr(
      Lambda4)));


   return twoLoop * beta_g2up;
}

/**
 * Calculates the 3-loop beta function of g2up.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2up_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2up;

   beta_g2up = 0;


   return threeLoop * beta_g2up;
}

/**
 * Calculates the 4-loop beta function of g2up.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2up_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2up;

   beta_g2up = 0;


   return fourLoop * beta_g2up;
}

/**
 * Calculates the 5-loop beta function of g2up.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_susy_parameters::calc_beta_g2up_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2up;

   beta_g2up = 0;


   return fiveLoop * beta_g2up;
}

} // namespace flexiblesusy
