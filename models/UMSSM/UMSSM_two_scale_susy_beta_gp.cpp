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

// File generated at Mon 27 Feb 2017 13:34:02

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of gp.
 *
 * @return one-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_gp_one_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);


   double beta_gp;

   beta_gp = Re(Power(gp,3)*oneOver16PiSqr*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr
      (QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr
      (Qv)));


   return beta_gp;
}

/**
 * Calculates the two-loop beta function of gp.
 *
 * @return two-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_gp_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_gp;

   beta_gp = Re(0.4*Power(gp,3)*twoLoop*(-10*AbsSqr(Lambdax)*(Sqr(QHd) +
      Sqr(QHu) + Sqr(Qs)) + 3*Sqr(g1)*(2*Sqr(Qd) + 6*Sqr(Qe) + Sqr(QHd) + Sqr(
      QHu) + 3*Sqr(Ql) + Sqr(Qq) + 8*Sqr(Qu)) + 5*(2*(9*Power(Qd,4) + 3*Power(
      Qe,4) + 2*Power(QHd,4) + 2*Power(QHu,4) + 6*Power(Ql,4) + 18*Power(Qq,4)
      + Power(Qs,4) + 9*Power(Qu,4) + 3*Power(Qv,4))*Sqr(gp) - 6*traceYdAdjYd*
      Sqr(Qd) - 2*traceYeAdjYe*Sqr(Qe) - 6*traceYdAdjYd*Sqr(QHd) - 2*
      traceYeAdjYe*Sqr(QHd) + 3*Sqr(g2)*Sqr(QHd) - 6*traceYuAdjYu*Sqr(QHu) - 2*
      traceYvAdjYv*Sqr(QHu) + 3*Sqr(g2)*Sqr(QHu) - 2*traceYeAdjYe*Sqr(Ql) - 2*
      traceYvAdjYv*Sqr(Ql) + 9*Sqr(g2)*Sqr(Ql) - 6*traceYdAdjYd*Sqr(Qq) - 6*
      traceYuAdjYu*Sqr(Qq) + 27*Sqr(g2)*Sqr(Qq) - 6*traceYuAdjYu*Sqr(Qu) + 24*
      Sqr(g3)*(Sqr(Qd) + 2*Sqr(Qq) + Sqr(Qu)) - 2*traceYvAdjYv*Sqr(Qv))));


   return beta_gp;
}

/**
 * Calculates the three-loop beta function of gp.
 *
 * @return three-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_gp_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_gp;

   beta_gp = 0;


   return beta_gp;
}

} // namespace flexiblesusy
