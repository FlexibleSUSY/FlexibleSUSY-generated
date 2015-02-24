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

// File generated at Tue 24 Feb 2015 17:35:14

#include "UMSSM_two_scale_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the one-loop beta function of Yu.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yu_one_loop(const Susy_traces& susy_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = oneOver16PiSqr*(Yu*(3*traceYuAdjYu + AbsSqr(Lambdax) -
      0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr(g3) - 2*
      Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(Qq) - 2*Sqr(gp)*Sqr(Qu)) + Yu*Yd.adjoint
      ()*Yd + 3*(Yu*Yu.adjoint()*Yu));


   return beta_Yu;
}

/**
 * Calculates the two-loop beta function of Yu.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yu_two_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = twoLoop*(Yu*(6.095555555555555*Power(g1,4) + 7.5*Power(g2,4)
      - 1.7777777777777777*Power(g3,4) + 8*Power(gp,4)*Power(QHu,4) + 40*Power
      (gp,4)*Power(Qq,4) + 22*Power(gp,4)*Power(Qu,4) - 3*traceYdAdjYuYuAdjYd -
      9*traceYuAdjYuYuAdjYu + Sqr(g1)*Sqr(g2) + 3.022222222222222*Sqr(g1)*Sqr(
      g3) + 8*Sqr(g2)*Sqr(g3) + 3.6*Qd*QHu*Sqr(g1)*Sqr(gp) + 3.6*Qe*QHu*Sqr(g1)
      *Sqr(gp) - 1.2*QHd*QHu*Sqr(g1)*Sqr(gp) - 3.6*QHu*Ql*Sqr(g1)*Sqr(gp) + 1.2
      *Qd*Qq*Sqr(g1)*Sqr(gp) + 1.2*Qe*Qq*Sqr(g1)*Sqr(gp) - 0.4*QHd*Qq*Sqr(g1)*
      Sqr(gp) + 4*QHu*Qq*Sqr(g1)*Sqr(gp) - 1.2*Ql*Qq*Sqr(g1)*Sqr(gp) - 4.8*Qd*
      Qu*Sqr(g1)*Sqr(gp) - 4.8*Qe*Qu*Sqr(g1)*Sqr(gp) + 1.6*QHd*Qu*Sqr(g1)*Sqr(
      gp) - 8.8*QHu*Qu*Sqr(g1)*Sqr(gp) + 4.8*Ql*Qu*Sqr(g1)*Sqr(gp) - 7.2*Qq*Qu*
      Sqr(g1)*Sqr(gp) + 2.4*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 6*Sqr(g2)*Sqr(gp)*Sqr(
      QHu) + 18*Power(gp,4)*Sqr(Qd)*Sqr(QHu) + 6*Power(gp,4)*Sqr(Qe)*Sqr(QHu) +
      4*Power(gp,4)*Sqr(QHd)*Sqr(QHu) + 12*Power(gp,4)*Sqr(QHu)*Sqr(Ql) +
      1.3333333333333333*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 6*Sqr(g2)*Sqr(gp)*Sqr(Qq) +
      10.666666666666666*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 18*Power(gp,4)*Sqr(Qd)*Sqr(
      Qq) + 6*Power(gp,4)*Sqr(Qe)*Sqr(Qq) + 4*Power(gp,4)*Sqr(QHd)*Sqr(Qq) + 40
      *Power(gp,4)*Sqr(QHu)*Sqr(Qq) + 12*Power(gp,4)*Sqr(Ql)*Sqr(Qq) + 2*Power(
      gp,4)*Sqr(QHu)*Sqr(Qs) + 2*Power(gp,4)*Sqr(Qq)*Sqr(Qs) + AbsSqr(Lambdax)*
      (-3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) + Sqr(Qs
      ))) + 11.733333333333333*Sqr(g1)*Sqr(gp)*Sqr(Qu) + 10.666666666666666*Sqr
      (g3)*Sqr(gp)*Sqr(Qu) + 18*Power(gp,4)*Sqr(Qd)*Sqr(Qu) + 6*Power(gp,4)*Sqr
      (Qe)*Sqr(Qu) + 4*Power(gp,4)*Sqr(QHd)*Sqr(Qu) + 22*Power(gp,4)*Sqr(QHu)*
      Sqr(Qu) + 12*Power(gp,4)*Sqr(Ql)*Sqr(Qu) + 54*Power(gp,4)*Sqr(Qq)*Sqr(Qu)
      + 2*Power(gp,4)*Sqr(Qs)*Sqr(Qu) + 0.4*traceYuAdjYu*(2*Sqr(g1) + 5*(8*Sqr
      (g3) + 3*Sqr(gp)*(-Sqr(QHu) + Sqr(Qq) + Sqr(Qu)))) - 3*Sqr(Conj(Lambdax))
      *Sqr(Lambdax)) + (-3*traceYdAdjYd - traceYeAdjYe - AbsSqr(Lambdax) + 0.4*
      Sqr(g1) + 2*Sqr(gp)*Sqr(Qd) + 2*Sqr(gp)*Sqr(QHd) - 2*Sqr(gp)*Sqr(Qq))*(Yu
      *Yd.adjoint()*Yd) - 9*traceYuAdjYu*(Yu*Yu.adjoint()*Yu) - 3*AbsSqr(
      Lambdax)*(Yu*Yu.adjoint()*Yu) + 0.4*Sqr(g1)*(Yu*Yu.adjoint()*Yu) + 6*Sqr(
      g2)*(Yu*Yu.adjoint()*Yu) + 6*Sqr(gp)*Sqr(QHu)*(Yu*Yu.adjoint()*Yu) + 2*
      Sqr(gp)*Sqr(Qq)*(Yu*Yu.adjoint()*Yu) - 2*Sqr(gp)*Sqr(Qu)*(Yu*Yu.adjoint()
      *Yu) - 2*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yu*Yd.adjoint()*Yd*
      Yu.adjoint()*Yu) - 4*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu));


   return beta_Yu;
}

} // namespace flexiblesusy
