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

// File generated at Tue 7 Jul 2015 12:46:27

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


   double beta_vu;

   beta_vu = Re(0.1*oneOver16PiSqr*vu*(-30*traceYuAdjYu - 10*AbsSqr(
      Lambdax) + 3*Sqr(g1) + 15*Sqr(g2) + 20*Sqr(gp)*Sqr(QHu)));


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
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_vu;

   beta_vu = Re(-0.005*twoLoop*vu*(207*Power(g1,4) + 275*Power(g2,4) +
      800*Power(gp,4)*Power(QHu,4) - 600*traceYdAdjYuYuAdjYd - 1800*
      traceYuAdjYuYuAdjYu + 90*Sqr(g1)*Sqr(g2) + 360*Qd*QHu*Sqr(g1)*Sqr(gp) +
      360*Qe*QHu*Sqr(g1)*Sqr(gp) - 120*QHd*QHu*Sqr(g1)*Sqr(gp) - 360*QHu*Ql*Sqr
      (g1)*Sqr(gp) + 360*QHu*Qq*Sqr(g1)*Sqr(gp) - 720*QHu*Qu*Sqr(g1)*Sqr(gp) +
      240*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 600*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 1800*Power(
      gp,4)*Sqr(Qd)*Sqr(QHu) + 600*Power(gp,4)*Sqr(Qe)*Sqr(QHu) + 400*Power(gp,
      4)*Sqr(QHd)*Sqr(QHu) + 1200*Power(gp,4)*Sqr(QHu)*Sqr(Ql) + 3600*Power(gp,
      4)*Sqr(QHu)*Sqr(Qq) + 200*Power(gp,4)*Sqr(QHu)*Sqr(Qs) + 20*AbsSqr(
      Lambdax)*(-30*traceYdAdjYd - 10*traceYeAdjYe + 3*Sqr(g1) + 15*Sqr(g2) +
      20*Sqr(gp)*Sqr(QHd) + 20*Sqr(gp)*Sqr(Qs)) + 1800*Power(gp,4)*Sqr(QHu)*Sqr
      (Qu) + 20*traceYuAdjYu*(17*Sqr(g1) + 45*Sqr(g2) + 160*Sqr(g3) + 60*Sqr(gp
      )*Sqr(Qq) + 60*Sqr(gp)*Sqr(Qu)) - 600*Sqr(Conj(Lambdax))*Sqr(Lambdax)));


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
