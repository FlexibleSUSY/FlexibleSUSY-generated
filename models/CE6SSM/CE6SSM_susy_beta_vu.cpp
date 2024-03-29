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


#include "CE6SSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of vu.
 *
 * @return 1-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_vu_1_loop(const Susy_traces& susy_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_vu;

   beta_vu = Re(0.1*vu*(-30*traceYuAdjYu - 10*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*
      Sqr(g2) + 2*Sqr(gN)));


   return oneLoop * beta_vu;
}

/**
 * Calculates the 2-loop beta function of vu.
 *
 * @return 2-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_vu_2_loop(const Susy_traces& susy_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_vu;

   beta_vu = Re(-0.005*vu*(-600*traceYdAdjYuYuAdjYd - 1800*traceYuAdjYuYuAdjYu
      - 600*traceKappaAdjKappa*AbsSqr(Lambdax) - 400*traceLambda12AdjLambda12*
      AbsSqr(Lambdax) - 600*traceYdAdjYd*AbsSqr(Lambdax) - 200*traceYeAdjYe*
      AbsSqr(Lambdax) + 297*Quad(g1) + 725*Quad(g2) + 192*Quad(gN) + 340*
      traceYuAdjYu*Sqr(g1) + 60*AbsSqr(Lambdax)*Sqr(g1) + 900*traceYuAdjYu*Sqr(
      g2) + 300*AbsSqr(Lambdax)*Sqr(g2) + 90*Sqr(g1)*Sqr(g2) + 3200*
      traceYuAdjYu*Sqr(g3) + 60*traceYuAdjYu*Sqr(gN) + 340*AbsSqr(Lambdax)*Sqr(
      gN) + 36*Sqr(g1)*Sqr(gN) + 60*Sqr(g2)*Sqr(gN) - 600*Sqr(Conj(Lambdax))*
      Sqr(Lambdax)));


   return twoLoop * beta_vu;
}

/**
 * Calculates the 3-loop beta function of vu.
 *
 * @return 3-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_vu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return threeLoop * beta_vu;
}

/**
 * Calculates the 4-loop beta function of vu.
 *
 * @return 4-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_vu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return fourLoop * beta_vu;
}

/**
 * Calculates the 5-loop beta function of vu.
 *
 * @return 5-loop beta function
 */
double CE6SSM_susy_parameters::calc_beta_vu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return fiveLoop * beta_vu;
}

} // namespace flexiblesusy
