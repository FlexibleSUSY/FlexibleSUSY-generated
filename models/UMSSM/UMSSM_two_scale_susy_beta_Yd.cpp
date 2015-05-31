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
 * Calculates the one-loop beta function of Yd.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yd_one_loop(const Susy_traces& susy_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto QHd = INPUT(QHd);
   const auto Qq = INPUT(Qq);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(Yd*(3*traceYdAdjYd + traceYeAdjYe + AbsSqr(
      Lambdax) - 0.4666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr
      (g3) - 2*Sqr(gp)*Sqr(Qd) - 2*Sqr(gp)*Sqr(QHd) - 2*Sqr(gp)*Sqr(Qq)) + 3*(
      Yd*Yd.adjoint()*Yd) + Yd*Yu.adjoint()*Yu)).real();


   return beta_Yd;
}

/**
 * Calculates the two-loop beta function of Yd.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yd_two_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(Yd*(3.188888888888889*Power(g1,4) + 7.5*Power(g2,4
      ) - 1.7777777777777777*Power(g3,4) + 22*Power(gp,4)*Power(Qd,4) + 8*Power
      (gp,4)*Power(QHd,4) + 40*Power(gp,4)*Power(Qq,4) - 9*traceYdAdjYdYdAdjYd
      - 3*traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe + 1.2*traceYeAdjYe*Sqr(g1
      ) + Sqr(g1)*Sqr(g2) + 0.8888888888888888*Sqr(g1)*Sqr(g3) + 8*Sqr(g2)*Sqr(
      g3) + 2.4*Qd*Qe*Sqr(g1)*Sqr(gp) - 4.4*Qd*QHd*Sqr(g1)*Sqr(gp) - 3.6*Qe*QHd
      *Sqr(g1)*Sqr(gp) + 0.8*Qd*QHu*Sqr(g1)*Sqr(gp) - 1.2*QHd*QHu*Sqr(g1)*Sqr(
      gp) - 2.4*Qd*Ql*Sqr(g1)*Sqr(gp) + 3.6*QHd*Ql*Sqr(g1)*Sqr(gp) + 3.6*Qd*Qq*
      Sqr(g1)*Sqr(gp) + 1.2*Qe*Qq*Sqr(g1)*Sqr(gp) - 4*QHd*Qq*Sqr(g1)*Sqr(gp) +
      0.4*QHu*Qq*Sqr(g1)*Sqr(gp) - 1.2*Ql*Qq*Sqr(g1)*Sqr(gp) - 4.8*Qd*Qu*Sqr(g1
      )*Sqr(gp) + 7.2*QHd*Qu*Sqr(g1)*Sqr(gp) - 2.4*Qq*Qu*Sqr(g1)*Sqr(gp) +
      2.933333333333333*Sqr(g1)*Sqr(gp)*Sqr(Qd) + 10.666666666666666*Sqr(g3)*
      Sqr(gp)*Sqr(Qd) + 2*traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 6*Power(gp,4)*Sqr(Qd)*
      Sqr(Qe) - 2*traceYeAdjYe*Sqr(gp)*Sqr(QHd) + 2.4*Sqr(g1)*Sqr(gp)*Sqr(QHd)
      + 6*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 22*Power(gp,4)*Sqr(Qd)*Sqr(QHd) + 6*Power(
      gp,4)*Sqr(Qe)*Sqr(QHd) + 4*Power(gp,4)*Sqr(Qd)*Sqr(QHu) + 4*Power(gp,4)*
      Sqr(QHd)*Sqr(QHu) + 2*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 12*Power(gp,4)*Sqr(
      Qd)*Sqr(Ql) + 12*Power(gp,4)*Sqr(QHd)*Sqr(Ql) + 1.3333333333333333*Sqr(g1
      )*Sqr(gp)*Sqr(Qq) + 6*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 10.666666666666666*Sqr(g3
      )*Sqr(gp)*Sqr(Qq) + 54*Power(gp,4)*Sqr(Qd)*Sqr(Qq) + 6*Power(gp,4)*Sqr(Qe
      )*Sqr(Qq) + 40*Power(gp,4)*Sqr(QHd)*Sqr(Qq) + 4*Power(gp,4)*Sqr(QHu)*Sqr(
      Qq) + 12*Power(gp,4)*Sqr(Ql)*Sqr(Qq) - 0.4*traceYdAdjYd*(Sqr(g1) - 5*(8*
      Sqr(g3) + 3*Sqr(gp)*(Sqr(Qd) - Sqr(QHd) + Sqr(Qq)))) + 2*Power(gp,4)*Sqr(
      Qd)*Sqr(Qs) + 2*Power(gp,4)*Sqr(QHd)*Sqr(Qs) + 2*Power(gp,4)*Sqr(Qq)*Sqr(
      Qs) + AbsSqr(Lambdax)*(-3*traceYuAdjYu + 2*Sqr(gp)*(-Sqr(QHd) + Sqr(QHu)
      + Sqr(Qs))) + 18*Power(gp,4)*Sqr(Qd)*Sqr(Qu) + 18*Power(gp,4)*Sqr(QHd)*
      Sqr(Qu) + 18*Power(gp,4)*Sqr(Qq)*Sqr(Qu) - 3*Sqr(Conj(Lambdax))*Sqr(
      Lambdax)) + (-9*traceYdAdjYd - 3*traceYeAdjYe - 3*AbsSqr(Lambdax) + 0.8*
      Sqr(g1) + 6*Sqr(g2) - 2*Sqr(gp)*Sqr(Qd) + 6*Sqr(gp)*Sqr(QHd) + 2*Sqr(gp)*
      Sqr(Qq))*(Yd*Yd.adjoint()*Yd) - 3*traceYuAdjYu*(Yd*Yu.adjoint()*Yu) -
      AbsSqr(Lambdax)*(Yd*Yu.adjoint()*Yu) + 0.8*Sqr(g1)*(Yd*Yu.adjoint()*Yu) +
      2*Sqr(gp)*Sqr(QHu)*(Yd*Yu.adjoint()*Yu) - 2*Sqr(gp)*Sqr(Qq)*(Yd*
      Yu.adjoint()*Yu) + 2*Sqr(gp)*Sqr(Qu)*(Yd*Yu.adjoint()*Yu) - 4*(Yd*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd)
      - 2*(Yd*Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

} // namespace flexiblesusy
