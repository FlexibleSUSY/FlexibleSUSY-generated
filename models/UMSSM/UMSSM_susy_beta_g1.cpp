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

// File generated at Tue 22 Jan 2019 17:28:25

#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of g1.
 *
 * @return 1-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_g1_1_loop(const Susy_traces& susy_traces) const
{


   double beta_g1;

   beta_g1 = Re(6.6*oneOver16PiSqr*Cube(g1));


   return beta_g1;
}

/**
 * Calculates the 2-loop beta function of g1.
 *
 * @return 2-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_g1_2_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_g1;

   beta_g1 = Re(0.04*twoLoop*Cube(g1)*(-70*traceYdAdjYd - 90*traceYeAdjYe - 130
      *traceYuAdjYu - 30*traceYvAdjYv - 30*AbsSqr(Lambdax) + 199*Sqr(g1) + 135*
      Sqr(g2) + 440*Sqr(g3) + 60*Sqr(gp)*Sqr(Qd) + 180*Sqr(gp)*Sqr(Qe) + 30*Sqr
      (gp)*Sqr(QHd) + 30*Sqr(gp)*Sqr(QHu) + 90*Sqr(gp)*Sqr(Ql) + 30*Sqr(gp)*Sqr
      (Qq) + 240*Sqr(gp)*Sqr(Qu)));


   return beta_g1;
}

/**
 * Calculates the 3-loop beta function of g1.
 *
 * @return 3-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_g1_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1;

   beta_g1 = 0;


   return beta_g1;
}

/**
 * Calculates the 4-loop beta function of g1.
 *
 * @return 4-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_g1_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1;

   beta_g1 = 0;


   return beta_g1;
}

/**
 * Calculates the 5-loop beta function of g1.
 *
 * @return 5-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_g1_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_g1;

   beta_g1 = 0;


   return beta_g1;
}

} // namespace flexiblesusy
