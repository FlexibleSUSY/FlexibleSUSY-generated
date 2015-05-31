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

// File generated at Sun 31 May 2015 12:30:41

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Ye.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_one_loop(const Susy_traces& susy_traces) const
{
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(Ye*(3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 1.8*Sqr(g1) - 3*Sqr(g2) - 2*Sqr(gp)*Sqr(Qe) - 2*Sqr(gp)*Sqr(
      QHd) - 2*Sqr(gp)*Sqr(Ql)) + 3*(Ye*Ye.adjoint()*Ye))).real();


   return beta_Ye;
}

/**
 * Calculates the two-loop beta function of Ye.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_two_loop(const Susy_traces& susy_traces) const
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(0.1*Ye*(135*Power(g1,4) + 75*Power(g2,4) + 100*
      Power(gp,4)*Power(Qe,4) + 80*Power(gp,4)*Power(QHd,4) + 160*Power(gp,4)*
      Power(Ql,4) - 90*traceYdAdjYdYdAdjYd - 30*traceYdAdjYuYuAdjYd - 30*
      traceYeAdjYeYeAdjYe + 12*traceYeAdjYe*Sqr(g1) + 18*Sqr(g1)*Sqr(g2) + 72*
      Qd*Qe*Sqr(g1)*Sqr(gp) - 36*Qd*QHd*Sqr(g1)*Sqr(gp) - 60*Qe*QHd*Sqr(g1)*Sqr
      (gp) + 24*Qe*QHu*Sqr(g1)*Sqr(gp) - 12*QHd*QHu*Sqr(g1)*Sqr(gp) - 36*Qd*Ql*
      Sqr(g1)*Sqr(gp) - 108*Qe*Ql*Sqr(g1)*Sqr(gp) + 48*QHd*Ql*Sqr(g1)*Sqr(gp) -
      12*QHu*Ql*Sqr(g1)*Sqr(gp) + 72*Qe*Qq*Sqr(g1)*Sqr(gp) - 36*QHd*Qq*Sqr(g1)
      *Sqr(gp) - 36*Ql*Qq*Sqr(g1)*Sqr(gp) - 144*Qe*Qu*Sqr(g1)*Sqr(gp) + 72*QHd*
      Qu*Sqr(g1)*Sqr(gp) + 72*Ql*Qu*Sqr(g1)*Sqr(gp) + 20*traceYeAdjYe*Sqr(gp)*
      Sqr(Qe) + 120*Sqr(g1)*Sqr(gp)*Sqr(Qe) + 180*Power(gp,4)*Sqr(Qd)*Sqr(Qe) -
      20*traceYeAdjYe*Sqr(gp)*Sqr(QHd) + 24*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 60*Sqr(
      g2)*Sqr(gp)*Sqr(QHd) + 180*Power(gp,4)*Sqr(Qd)*Sqr(QHd) + 100*Power(gp,4)
      *Sqr(Qe)*Sqr(QHd) + 40*Power(gp,4)*Sqr(Qe)*Sqr(QHu) + 40*Power(gp,4)*Sqr(
      QHd)*Sqr(QHu) + 20*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 48*Sqr(g1)*Sqr(gp)*Sqr(
      Ql) + 60*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 180*Power(gp,4)*Sqr(Qd)*Sqr(Ql) + 180*
      Power(gp,4)*Sqr(Qe)*Sqr(Ql) + 160*Power(gp,4)*Sqr(QHd)*Sqr(Ql) + 40*Power
      (gp,4)*Sqr(QHu)*Sqr(Ql) + 360*Power(gp,4)*Sqr(Qe)*Sqr(Qq) + 360*Power(gp,
      4)*Sqr(QHd)*Sqr(Qq) + 360*Power(gp,4)*Sqr(Ql)*Sqr(Qq) - 4*traceYdAdjYd*(
      Sqr(g1) - 5*(8*Sqr(g3) + 3*Sqr(gp)*(Sqr(Qd) - Sqr(QHd) + Sqr(Qq)))) - 10*
      AbsSqr(Lambdax)*(3*traceYuAdjYu + 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) - Sqr(Qs
      ))) + 20*Power(gp,4)*Sqr(Qe)*Sqr(Qs) + 20*Power(gp,4)*Sqr(QHd)*Sqr(Qs) +
      20*Power(gp,4)*Sqr(Ql)*Sqr(Qs) + 180*Power(gp,4)*Sqr(Qe)*Sqr(Qu) + 180*
      Power(gp,4)*Sqr(QHd)*Sqr(Qu) + 180*Power(gp,4)*Sqr(Ql)*Sqr(Qu) - 30*Sqr(
      Conj(Lambdax))*Sqr(Lambdax)) + (-9*traceYdAdjYd - 3*traceYeAdjYe - 3*
      AbsSqr(Lambdax) + 6*Sqr(g2) - 2*Sqr(gp)*Sqr(Qe) + 6*Sqr(gp)*Sqr(QHd) + 2*
      Sqr(gp)*Sqr(Ql))*(Ye*Ye.adjoint()*Ye) - 4*(Ye*Ye.adjoint()*Ye*Ye.adjoint(
      )*Ye))).real();


   return beta_Ye;
}

} // namespace flexiblesusy
