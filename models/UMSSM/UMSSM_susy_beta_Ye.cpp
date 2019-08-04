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

// File generated at Sun 4 Aug 2019 19:33:34

#include "UMSSM_susy_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT susy_traces

/**
 * Calculates the 1-loop beta function of Ye.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_1_loop(const Susy_traces& susy_traces) const
{
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto Ql = INPUT(Ql);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(-0.2*Ye*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*
      AbsSqr(Lambdax) + 9*Sqr(g1) + 15*Sqr(g2) + 10*Sqr(gp)*Sqr(Qe) + 10*Sqr(gp
      )*Sqr(QHd) + 10*Sqr(gp)*Sqr(Ql)) + 3*(Ye*Ye.adjoint()*Ye) + Ye*Yv.
      conjugate()*Yv.transpose())).real();


   return beta_Ye;
}

/**
 * Calculates the 2-loop beta function of Ye.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_2_loop(const Susy_traces& susy_traces) const
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
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(0.1*Ye*(-90*traceYdAdjYdYdAdjYd - 30*traceYdAdjYuYuAdjYd
       - 30*traceYeAdjYeYeAdjYe - 10*traceYvAdjYvTpYeconjYe - 30*traceYuAdjYu*
      AbsSqr(Lambdax) - 10*traceYvAdjYv*AbsSqr(Lambdax) + 135*Quad(g1) + 75*
      Quad(g2) + 100*Quad(gp)*Quad(Qe) + 80*Quad(gp)*Quad(QHd) + 160*Quad(gp)*
      Quad(Ql) - 4*traceYdAdjYd*Sqr(g1) + 12*traceYeAdjYe*Sqr(g1) + 18*Sqr(g1)*
      Sqr(g2) + 160*traceYdAdjYd*Sqr(g3) + 72*Qd*Qe*Sqr(g1)*Sqr(gp) - 36*Qd*QHd
      *Sqr(g1)*Sqr(gp) - 60*Qe*QHd*Sqr(g1)*Sqr(gp) + 24*Qe*QHu*Sqr(g1)*Sqr(gp)
      - 12*QHd*QHu*Sqr(g1)*Sqr(gp) - 36*Qd*Ql*Sqr(g1)*Sqr(gp) - 108*Qe*Ql*Sqr(
      g1)*Sqr(gp) + 48*QHd*Ql*Sqr(g1)*Sqr(gp) - 12*QHu*Ql*Sqr(g1)*Sqr(gp) + 72*
      Qe*Qq*Sqr(g1)*Sqr(gp) - 36*QHd*Qq*Sqr(g1)*Sqr(gp) - 36*Ql*Qq*Sqr(g1)*Sqr(
      gp) - 144*Qe*Qu*Sqr(g1)*Sqr(gp) + 72*QHd*Qu*Sqr(g1)*Sqr(gp) + 72*Ql*Qu*
      Sqr(g1)*Sqr(gp) + 60*traceYdAdjYd*Sqr(gp)*Sqr(Qd) + 20*traceYeAdjYe*Sqr(
      gp)*Sqr(Qe) + 120*Sqr(g1)*Sqr(gp)*Sqr(Qe) + 180*Quad(gp)*Sqr(Qd)*Sqr(Qe)
      - 60*traceYdAdjYd*Sqr(gp)*Sqr(QHd) - 20*traceYeAdjYe*Sqr(gp)*Sqr(QHd) -
      20*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 24*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 60*
      Sqr(g2)*Sqr(gp)*Sqr(QHd) + 180*Quad(gp)*Sqr(Qd)*Sqr(QHd) + 100*Quad(gp)*
      Sqr(Qe)*Sqr(QHd) + 20*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 40*Quad(gp)*Sqr(
      Qe)*Sqr(QHu) + 40*Quad(gp)*Sqr(QHd)*Sqr(QHu) + 20*traceYeAdjYe*Sqr(gp)*
      Sqr(Ql) + 48*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 60*Sqr(g2)*Sqr(gp)*Sqr(Ql) + 180*
      Quad(gp)*Sqr(Qd)*Sqr(Ql) + 180*Quad(gp)*Sqr(Qe)*Sqr(Ql) + 160*Quad(gp)*
      Sqr(QHd)*Sqr(Ql) + 40*Quad(gp)*Sqr(QHu)*Sqr(Ql) + 60*traceYdAdjYd*Sqr(gp)
      *Sqr(Qq) + 360*Quad(gp)*Sqr(Qe)*Sqr(Qq) + 360*Quad(gp)*Sqr(QHd)*Sqr(Qq) +
      360*Quad(gp)*Sqr(Ql)*Sqr(Qq) + 20*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*
      Quad(gp)*Sqr(Qe)*Sqr(Qs) + 20*Quad(gp)*Sqr(QHd)*Sqr(Qs) + 20*Quad(gp)*Sqr
      (Ql)*Sqr(Qs) + 180*Quad(gp)*Sqr(Qe)*Sqr(Qu) + 180*Quad(gp)*Sqr(QHd)*Sqr(
      Qu) + 180*Quad(gp)*Sqr(Ql)*Sqr(Qu) + 60*Quad(gp)*Sqr(Qe)*Sqr(Qv) + 60*
      Quad(gp)*Sqr(QHd)*Sqr(Qv) + 60*Quad(gp)*Sqr(Ql)*Sqr(Qv) - 30*Sqr(Conj(
      Lambdax))*Sqr(Lambdax)) + (-9*traceYdAdjYd - 3*traceYeAdjYe - 3*AbsSqr(
      Lambdax) + 6*Sqr(g2) - 2*Sqr(gp)*Sqr(Qe) + 6*Sqr(gp)*Sqr(QHd) + 2*Sqr(gp)
      *Sqr(Ql))*(Ye*Ye.adjoint()*Ye) + (-3*traceYuAdjYu - traceYvAdjYv - AbsSqr
      (Lambdax) + 2*Sqr(gp)*Sqr(QHu) - 2*Sqr(gp)*Sqr(Ql) + 2*Sqr(gp)*Sqr(Qv))*(
      Ye*Yv.conjugate()*Yv.transpose()) - 4*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye
      ) - 2*(Ye*Yv.conjugate()*Yv.transpose()*Ye.adjoint()*Ye) - 2*(Ye*Yv.
      conjugate()*Yv.transpose()*Yv.conjugate()*Yv.transpose()))).real();


   return beta_Ye;
}

/**
 * Calculates the 3-loop beta function of Ye.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

/**
 * Calculates the 4-loop beta function of Ye.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

/**
 * Calculates the 5-loop beta function of Ye.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Ye_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = ZEROMATRIX(3,3);


   return beta_Ye;
}

} // namespace flexiblesusy
