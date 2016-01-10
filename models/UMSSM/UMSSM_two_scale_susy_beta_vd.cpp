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

// File generated at Sun 10 Jan 2016 15:33:43

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

   beta_vd = Re(0.1*oneOver16PiSqr*vd*(-30*traceYdAdjYd - 10*traceYeAdjYe
      - 10*AbsSqr(Lambdax) + 3*Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gp)*Sqr(QHd)));


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

   beta_vd = Re(-0.005*twoLoop*vd*(261*Power(g1,4) + 275*Power(g2,4) +
      800*Power(gp,4)*Power(QHd,4) - 1800*traceYdAdjYdYdAdjYd - 600*
      traceYdAdjYuYuAdjYd - 600*traceYeAdjYeYeAdjYe - 200*
      traceYvAdjYvTpYeconjYe + 300*traceYeAdjYe*Sqr(g1) + 300*traceYeAdjYe*Sqr(
      g2) + 90*Sqr(g1)*Sqr(g2) - 360*Qd*QHd*Sqr(g1)*Sqr(gp) - 360*Qe*QHd*Sqr(g1
      )*Sqr(gp) - 120*QHd*QHu*Sqr(g1)*Sqr(gp) + 360*QHd*Ql*Sqr(g1)*Sqr(gp) -
      360*QHd*Qq*Sqr(g1)*Sqr(gp) + 720*QHd*Qu*Sqr(g1)*Sqr(gp) - 360*QHd*Qv*Sqr(
      g1)*Sqr(gp) + 400*traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 240*Sqr(g1)*Sqr(gp)*Sqr(
      QHd) + 600*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 1800*Power(gp,4)*Sqr(Qd)*Sqr(QHd) +
      600*Power(gp,4)*Sqr(Qe)*Sqr(QHd) + 400*Power(gp,4)*Sqr(QHd)*Sqr(QHu) +
      400*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 1200*Power(gp,4)*Sqr(QHd)*Sqr(Ql) +
      3600*Power(gp,4)*Sqr(QHd)*Sqr(Qq) + 100*traceYdAdjYd*(Sqr(g1) + 9*Sqr(g2)
      + 32*Sqr(g3) + 12*Sqr(gp)*Sqr(Qd) + 12*Sqr(gp)*Sqr(Qq)) + 200*Power(gp,4
      )*Sqr(QHd)*Sqr(Qs) + 20*AbsSqr(Lambdax)*(-30*traceYuAdjYu - 10*
      traceYvAdjYv + 3*Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu) + 20*Sqr(gp)*
      Sqr(Qs)) + 1800*Power(gp,4)*Sqr(QHd)*Sqr(Qu) + 600*Power(gp,4)*Sqr(QHd)*
      Sqr(Qv) - 600*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


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
