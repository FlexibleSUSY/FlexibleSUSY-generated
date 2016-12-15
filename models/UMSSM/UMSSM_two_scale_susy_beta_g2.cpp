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

// File generated at Thu 15 Dec 2016 12:51:26

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of g2.
 *
 * @return one-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_g2_one_loop(const Susy_traces& susy_traces) const
{


   double beta_g2;

   beta_g2 = Re(Power(g2,3)*oneOver16PiSqr);


   return beta_g2;
}

/**
 * Calculates the two-loop beta function of g2.
 *
 * @return two-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_g2_two_loop(const Susy_traces& susy_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_g2;

   beta_g2 = Re(0.2*Power(g2,3)*twoLoop*(-10*AbsSqr(Lambdax) + 9*Sqr(g1)
      + 5*(-6*traceYdAdjYd - 2*traceYeAdjYe - 6*traceYuAdjYu - 2*traceYvAdjYv +
      25*Sqr(g2) + 24*Sqr(g3) + 2*Sqr(gp)*Sqr(QHd) + 2*Sqr(gp)*Sqr(QHu) + 6*
      Sqr(gp)*Sqr(Ql) + 18*Sqr(gp)*Sqr(Qq))));


   return beta_g2;
}

/**
 * Calculates the three-loop beta function of g2.
 *
 * @return three-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_g2_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g2;

   beta_g2 = 0;


   return beta_g2;
}

} // namespace flexiblesusy
