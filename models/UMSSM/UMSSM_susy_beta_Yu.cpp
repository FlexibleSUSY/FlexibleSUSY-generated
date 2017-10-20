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

// File generated at Fri 20 Oct 2017 08:51:42

#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yu.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yu_1_loop(const Susy_traces& susy_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(Yu*(3*traceYuAdjYu + traceYvAdjYv + AbsSqr(
      Lambdax) - 0.8666666666666667*Sqr(g1) - 3*Sqr(g2) - 5.333333333333333*Sqr
      (g3) - 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(Qq) - 2*Sqr(gp)*Sqr(Qu)) + Yu*
      Yd.adjoint()*Yd + 3*(Yu*Yu.adjoint()*Yu))).real();


   return beta_Yu;
}

/**
 * Calculates the 2-loop beta function of Yu.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yu_2_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0022222222222222222*Yu*(2743*Quad(g1) + 450*
      AbsSqr(Lambdax)*(-3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(QHd) -
      Sqr(QHu) + Sqr(Qs))) + 10*Sqr(g1)*(45*Sqr(g2) + 2*(68*Sqr(g3) + 3*(6*
      traceYuAdjYu + Sqr(gp)*(-9*QHd*QHu - 27*QHu*Ql - 3*QHd*Qq + 30*QHu*Qq - 9
      *Ql*Qq + 9*Qd*(3*QHu + Qq - 4*Qu) + 9*Qe*(3*QHu + Qq - 4*Qu) + 12*QHd*Qu
      - 66*QHu*Qu + 36*Ql*Qu - 54*Qq*Qu + 18*Sqr(QHu) + 10*Sqr(Qq) + 88*Sqr(Qu)
      )))) + 25*(135*Quad(g2) + 36*Sqr(g2)*(4*Sqr(g3) + 3*Sqr(gp)*(Sqr(QHu) +
      Sqr(Qq))) + 2*(-16*Quad(g3) + 48*Sqr(g3)*(3*traceYuAdjYu + 2*Sqr(gp)*(Sqr
      (Qq) + Sqr(Qu))) + 9*(-3*traceYdAdjYuYuAdjYd - 9*traceYuAdjYuYuAdjYu -
      traceYvAdjYvTpYeconjYe - 3*traceYvAdjYvYvAdjYv + 2*Sqr(gp)*(-((3*
      traceYuAdjYu + traceYvAdjYv)*Sqr(QHu)) + traceYvAdjYv*Sqr(Ql) + 3*
      traceYuAdjYu*Sqr(Qq) + 3*traceYuAdjYu*Sqr(Qu) + traceYvAdjYv*Sqr(Qv)) + 2
      *Quad(gp)*(4*Quad(QHu) + 20*Quad(Qq) + 11*Quad(Qu) + 2*Sqr(QHd)*Sqr(QHu)
      + 6*Sqr(QHu)*Sqr(Ql) + 2*Sqr(QHd)*Sqr(Qq) + 20*Sqr(QHu)*Sqr(Qq) + 6*Sqr(
      Ql)*Sqr(Qq) + Sqr(QHu)*Sqr(Qs) + Sqr(Qq)*Sqr(Qs) + 2*Sqr(QHd)*Sqr(Qu) +
      11*Sqr(QHu)*Sqr(Qu) + 6*Sqr(Ql)*Sqr(Qu) + 27*Sqr(Qq)*Sqr(Qu) + Sqr(Qs)*
      Sqr(Qu) + 9*Sqr(Qd)*(Sqr(QHu) + Sqr(Qq) + Sqr(Qu)) + 3*Sqr(Qe)*(Sqr(QHu)
      + Sqr(Qq) + Sqr(Qu)) + 3*Sqr(QHu)*Sqr(Qv) + 3*Sqr(Qq)*Sqr(Qv) + 3*Sqr(Qu)
      *Sqr(Qv))))) - 1350*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + (-3*traceYdAdjYd -
      traceYeAdjYe - AbsSqr(Lambdax) + 0.4*Sqr(g1) + 2*Sqr(gp)*Sqr(Qd) + 2*Sqr
      (gp)*Sqr(QHd) - 2*Sqr(gp)*Sqr(Qq))*(Yu*Yd.adjoint()*Yd) + (-9*
      traceYuAdjYu - 3*traceYvAdjYv - 3*AbsSqr(Lambdax) + 0.4*Sqr(g1) + 6*Sqr(
      g2) + 6*Sqr(gp)*Sqr(QHu) + 2*Sqr(gp)*Sqr(Qq) - 2*Sqr(gp)*Sqr(Qu))*(Yu*
      Yu.adjoint()*Yu) - 2*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yu*
      Yd.adjoint()*Yd*Yu.adjoint()*Yu) - 4*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu)
      )).real();


   return beta_Yu;
}

/**
 * Calculates the 3-loop beta function of Yu.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yu_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return beta_Yu;
}

} // namespace flexiblesusy
