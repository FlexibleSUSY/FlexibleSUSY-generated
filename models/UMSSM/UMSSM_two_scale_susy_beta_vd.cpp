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

// File generated at Tue 5 Sep 2017 11:43:46

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of vd.
 *
 * @return one-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vd_one_loop(const Susy_traces& susy_traces) const
{
   const auto QHd = INPUT(QHd);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_vd;

   beta_vd = Re(0.1*oneOver16PiSqr*vd*(-10*AbsSqr(Lambdax) + 3*Sqr(g1) +
      5*(-6*traceYdAdjYd - 2*traceYeAdjYe + 3*Sqr(g2) + 4*Sqr(gp)*Sqr(QHd))));


   return beta_vd;
}

/**
 * Calculates the two-loop beta function of vd.
 *
 * @return two-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vd_two_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qe = INPUT(Qe);
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   double beta_vd;

   beta_vd = Re(-0.005*twoLoop*vd*(207*Power(g1,4) + 10*Sqr(g1)*(10*(
      traceYdAdjYd + 3*traceYeAdjYe) + 9*Sqr(g2) - 12*QHd*(3*Qd + 3*Qe - 2*QHd
      + QHu - 3*Ql + 3*Qq - 6*Qu)*Sqr(gp)) + 20*AbsSqr(Lambdax)*(3*Sqr(g1) + 5*
      (-2*(3*traceYuAdjYu + traceYvAdjYv) + 3*Sqr(g2) + 4*Sqr(gp)*(Sqr(QHu) +
      Sqr(Qs)))) + 25*(11*Power(g2,4) + 12*Sqr(g2)*(3*traceYdAdjYd +
      traceYeAdjYe + 2*Sqr(gp)*Sqr(QHd)) + 8*(-9*traceYdAdjYdYdAdjYd - 3*
      traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe - traceYvAdjYvTpYeconjYe + 16
      *traceYdAdjYd*Sqr(g3) + 2*Sqr(gp)*(3*traceYdAdjYd*Sqr(Qd) + traceYeAdjYe*
      (Sqr(Qe) + Sqr(Ql)) + 3*traceYdAdjYd*Sqr(Qq)) + Power(gp,4)*Sqr(QHd)*(9*
      Sqr(Qd) + 3*Sqr(Qe) + 4*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) +
      Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv)))) - 600*Sqr(Conj(Lambdax))*Sqr(Lambdax))
      );


   return beta_vd;
}

/**
 * Calculates the three-loop beta function of vd.
 *
 * @return three-loop beta function
 */
double UMSSM_susy_parameters::calc_beta_vd_three_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_vd;

   beta_vd = 0;


   return beta_vd;
}

} // namespace flexiblesusy
