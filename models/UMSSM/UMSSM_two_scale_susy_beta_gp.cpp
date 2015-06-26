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

// File generated at Fri 26 Jun 2015 19:03:24

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


   double beta_gp;

   beta_gp = Re(Power(gp,3)*oneOver16PiSqr*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr
      (QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu)));


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
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_gp;

   beta_gp = Re(0.4*Power(gp,3)*twoLoop*(90*Power(Qd,4)*Sqr(gp) + 30*
      Power(Qe,4)*Sqr(gp) + 20*Power(QHd,4)*Sqr(gp) + 20*Power(QHu,4)*Sqr(gp) +
      60*Power(Ql,4)*Sqr(gp) + 180*Power(Qq,4)*Sqr(gp) + 10*Power(Qs,4)*Sqr(gp
      ) + 90*Power(Qu,4)*Sqr(gp) + 6*Sqr(g1)*Sqr(Qd) + 120*Sqr(g3)*Sqr(Qd) - 10
      *traceYeAdjYe*Sqr(Qe) + 18*Sqr(g1)*Sqr(Qe) - 10*traceYeAdjYe*Sqr(QHd) + 3
      *Sqr(g1)*Sqr(QHd) + 15*Sqr(g2)*Sqr(QHd) - 30*traceYuAdjYu*Sqr(QHu) + 3*
      Sqr(g1)*Sqr(QHu) + 15*Sqr(g2)*Sqr(QHu) - 10*traceYeAdjYe*Sqr(Ql) + 9*Sqr(
      g1)*Sqr(Ql) + 45*Sqr(g2)*Sqr(Ql) - 30*traceYuAdjYu*Sqr(Qq) + 3*Sqr(g1)*
      Sqr(Qq) + 135*Sqr(g2)*Sqr(Qq) + 240*Sqr(g3)*Sqr(Qq) - 30*traceYdAdjYd*(
      Sqr(Qd) + Sqr(QHd) + Sqr(Qq)) - 10*AbsSqr(Lambdax)*(Sqr(QHd) + Sqr(QHu) +
      Sqr(Qs)) - 30*traceYuAdjYu*Sqr(Qu) + 24*Sqr(g1)*Sqr(Qu) + 120*Sqr(g3)*
      Sqr(Qu)));


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
