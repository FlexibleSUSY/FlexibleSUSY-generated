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

// File generated at Tue 10 Oct 2017 22:13:20

#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Yd.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yd_1_loop(const Susy_traces& susy_traces) const
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
 * Calculates the 2-loop beta function of Yd.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yd_2_loop(const Susy_traces& susy_traces) const
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(0.011111111111111112*Yd*(287*Quad(g1) + 2*Sqr(g1)*
      (45*Sqr(g2) + 2*(20*Sqr(g3) + 3*(-3*traceYdAdjYd + 9*traceYeAdjYe + Sqr(
      gp)*(-9*QHd*QHu + 27*QHd*Ql - 30*QHd*Qq + 3*QHu*Qq - 9*Ql*Qq + 9*Qe*(-3*
      QHd + Qq) + 3*Qd*(6*Qe - 11*QHd + 2*QHu - 6*Ql + 9*Qq - 12*Qu) + 54*QHd*
      Qu - 18*Qq*Qu + 22*Sqr(Qd) + 18*Sqr(QHd) + 10*Sqr(Qq))))) - 90*AbsSqr(
      Lambdax)*(3*traceYuAdjYu + traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHd) - Sqr(QHu)
      - Sqr(Qs))) + 5*(135*Quad(g2) + 36*Sqr(g2)*(4*Sqr(g3) + 3*Sqr(gp)*(Sqr(
      QHd) + Sqr(Qq))) + 2*(-16*Quad(g3) + 48*Sqr(g3)*(3*traceYdAdjYd + 2*Sqr(
      gp)*(Sqr(Qd) + Sqr(Qq))) + 9*(-9*traceYdAdjYdYdAdjYd - 3*
      traceYdAdjYuYuAdjYd - 3*traceYeAdjYeYeAdjYe - traceYvAdjYvTpYeconjYe + 2*
      Sqr(gp)*(3*traceYdAdjYd*Sqr(Qd) + traceYeAdjYe*Sqr(Qe) - (3*traceYdAdjYd
      + traceYeAdjYe)*Sqr(QHd) + traceYeAdjYe*Sqr(Ql) + 3*traceYdAdjYd*Sqr(Qq))
      + 2*Quad(gp)*(11*Quad(Qd) + 4*Quad(QHd) + 20*Quad(Qq) + 2*Sqr(QHd)*Sqr(
      QHu) + 6*Sqr(QHd)*Sqr(Ql) + 20*Sqr(QHd)*Sqr(Qq) + 2*Sqr(QHu)*Sqr(Qq) + 6*
      Sqr(Ql)*Sqr(Qq) + 3*Sqr(Qe)*(Sqr(QHd) + Sqr(Qq)) + Sqr(QHd)*Sqr(Qs) + Sqr
      (Qq)*Sqr(Qs) + 9*Sqr(QHd)*Sqr(Qu) + 9*Sqr(Qq)*Sqr(Qu) + 3*Sqr(QHd)*Sqr(Qv
      ) + 3*Sqr(Qq)*Sqr(Qv) + Sqr(Qd)*(3*Sqr(Qe) + 11*Sqr(QHd) + 2*Sqr(QHu) + 6
      *Sqr(Ql) + 27*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv)))))) - 270*Sqr(
      Conj(Lambdax))*Sqr(Lambdax)) + (-9*traceYdAdjYd - 3*traceYeAdjYe - 3*
      AbsSqr(Lambdax) + 0.8*Sqr(g1) + 6*Sqr(g2) - 2*Sqr(gp)*Sqr(Qd) + 6*Sqr(gp)
      *Sqr(QHd) + 2*Sqr(gp)*Sqr(Qq))*(Yd*Yd.adjoint()*Yd) + (-3*traceYuAdjYu -
      traceYvAdjYv - AbsSqr(Lambdax) + 0.8*Sqr(g1) + 2*Sqr(gp)*Sqr(QHu) - 2*Sqr
      (gp)*Sqr(Qq) + 2*Sqr(gp)*Sqr(Qu))*(Yd*Yu.adjoint()*Yu) - 4*(Yd*Yd.adjoint
      ()*Yd*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(Yd*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();


   return beta_Yd;
}

/**
 * Calculates the 3-loop beta function of Yd.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yd_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy
