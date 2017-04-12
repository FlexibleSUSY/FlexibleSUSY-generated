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

// File generated at Wed 12 Apr 2017 12:25:18

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vu.
 *
 * @return one-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vu_one_loop(const Susy_traces& susy_traces) const
{
   const auto QHu = INPUT(QHu);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   double beta_vu;

   beta_vu = Re(0.1*oneOver16PiSqr*vu*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) +
      5*(-6*traceYuAdjYu - 2*traceYvAdjYv + 3*Sqr(g2) + 4*Sqr(gp)*Sqr(QHu))));


   return beta_vu;
}

/**
 * Calculates the two-loop beta function of vu.
 *
 * @return two-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vu_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHu = INPUT(QHu);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   double beta_vu;

   beta_vu = Re(-0.005*twoLoop*vu*(207*Power(g1,4) + 10*Sqr(g1)*(34*
      traceYuAdjYu + 6*traceYvAdjYv + 9*Sqr(g2) + 12*QHu*(3*Qd + 3*Qe - QHd + 2
      *QHu - 3*Ql + 3*Qq - 6*Qu)*Sqr(gp)) + 20*AbsSqr(Lambdax)*(3*Sqr(g1) + 5*(
      -2*(3*traceYdAdjYd + traceYeAdjYe) + 3*Sqr(g2) + 4*Sqr(gp)*(Sqr(QHd) +
      Sqr(Qs)))) + 25*(11*Power(g2,4) + 12*Sqr(g2)*(3*traceYuAdjYu +
      traceYvAdjYv + 2*Sqr(gp)*Sqr(QHu)) + 8*(-3*traceYdAdjYuYuAdjYd - 9*
      traceYuAdjYuYuAdjYu - traceYvAdjYvTpYeconjYe - 3*traceYvAdjYvYvAdjYv + 16
      *traceYuAdjYu*Sqr(g3) + Power(gp,4)*Sqr(QHu)*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*
      Sqr(QHd) + 4*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*
      Sqr(Qv)) + 2*Sqr(gp)*(3*traceYuAdjYu*Sqr(Qq) + 3*traceYuAdjYu*Sqr(Qu) +
      traceYvAdjYv*(Sqr(Ql) + Sqr(Qv))))) - 600*Sqr(Conj(Lambdax))*Sqr(Lambdax)
      ));


   return beta_vu;
}

/**
 * Calculates the three-loop beta function of vu.
 *
 * @return three-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vu_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vu;

   beta_vu = 0;


   return beta_vu;
}

} // namespace flexiblesusy
