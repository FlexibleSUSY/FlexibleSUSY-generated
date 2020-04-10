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

// File generated at Fri 10 Apr 2020 20:17:57

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

   beta_Yd = (oneOver16PiSqr*(-0.06666666666666667*Yd*(-45*traceYdAdjYd - 15*
      traceYeAdjYe - 15*AbsSqr(Lambdax) + 7*Sqr(g1) + 45*Sqr(g2) + 80*Sqr(g3) +
      30*Sqr(gp)*Sqr(Qd) + 30*Sqr(gp)*Sqr(QHd) + 30*Sqr(gp)*Sqr(Qq)) + 3*(Yd*Yd
      .adjoint()*Yd) + Yd*Yu.adjoint()*Yu)).real();


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
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(0.011111111111111112*Yd*(-810*traceYdAdjYdYdAdjYd - 270*
      traceYdAdjYuYuAdjYd - 270*traceYeAdjYeYeAdjYe - 90*traceYvAdjYvTpYeconjYe
       - 270*traceYuAdjYu*AbsSqr(Lambdax) - 90*traceYvAdjYv*AbsSqr(Lambdax) +
      287*Quad(g1) + 675*Quad(g2) - 160*Quad(g3) + 1980*Quad(gp)*Quad(Qd) + 720
      *Quad(gp)*Quad(QHd) + 3600*Quad(gp)*Quad(Qq) - 36*traceYdAdjYd*Sqr(g1) +
      108*traceYeAdjYe*Sqr(g1) + 90*Sqr(g1)*Sqr(g2) + 1440*traceYdAdjYd*Sqr(g3)
      + 80*Sqr(g1)*Sqr(g3) + 720*Sqr(g2)*Sqr(g3) + 216*Qd*Qe*Sqr(g1)*Sqr(gp) -
      396*Qd*QHd*Sqr(g1)*Sqr(gp) - 324*Qe*QHd*Sqr(g1)*Sqr(gp) + 72*Qd*QHu*Sqr(
      g1)*Sqr(gp) - 108*QHd*QHu*Sqr(g1)*Sqr(gp) - 216*Qd*Ql*Sqr(g1)*Sqr(gp) +
      324*QHd*Ql*Sqr(g1)*Sqr(gp) + 324*Qd*Qq*Sqr(g1)*Sqr(gp) + 108*Qe*Qq*Sqr(g1
      )*Sqr(gp) - 360*QHd*Qq*Sqr(g1)*Sqr(gp) + 36*QHu*Qq*Sqr(g1)*Sqr(gp) - 108*
      Ql*Qq*Sqr(g1)*Sqr(gp) - 432*Qd*Qu*Sqr(g1)*Sqr(gp) + 648*QHd*Qu*Sqr(g1)*
      Sqr(gp) - 216*Qq*Qu*Sqr(g1)*Sqr(gp) + 540*traceYdAdjYd*Sqr(gp)*Sqr(Qd) +
      264*Sqr(g1)*Sqr(gp)*Sqr(Qd) + 960*Sqr(g3)*Sqr(gp)*Sqr(Qd) + 180*
      traceYeAdjYe*Sqr(gp)*Sqr(Qe) + 540*Quad(gp)*Sqr(Qd)*Sqr(Qe) - 540*
      traceYdAdjYd*Sqr(gp)*Sqr(QHd) - 180*traceYeAdjYe*Sqr(gp)*Sqr(QHd) - 180*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 216*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 540*Sqr
      (g2)*Sqr(gp)*Sqr(QHd) + 1980*Quad(gp)*Sqr(Qd)*Sqr(QHd) + 540*Quad(gp)*Sqr
      (Qe)*Sqr(QHd) + 180*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 360*Quad(gp)*Sqr(
      Qd)*Sqr(QHu) + 360*Quad(gp)*Sqr(QHd)*Sqr(QHu) + 180*traceYeAdjYe*Sqr(gp)*
      Sqr(Ql) + 1080*Quad(gp)*Sqr(Qd)*Sqr(Ql) + 1080*Quad(gp)*Sqr(QHd)*Sqr(Ql)
      + 540*traceYdAdjYd*Sqr(gp)*Sqr(Qq) + 120*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 540*
      Sqr(g2)*Sqr(gp)*Sqr(Qq) + 960*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 4860*Quad(gp)*Sqr
      (Qd)*Sqr(Qq) + 540*Quad(gp)*Sqr(Qe)*Sqr(Qq) + 3600*Quad(gp)*Sqr(QHd)*Sqr(
      Qq) + 360*Quad(gp)*Sqr(QHu)*Sqr(Qq) + 1080*Quad(gp)*Sqr(Ql)*Sqr(Qq) + 180
      *AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 180*Quad(gp)*Sqr(Qd)*Sqr(Qs) + 180*
      Quad(gp)*Sqr(QHd)*Sqr(Qs) + 180*Quad(gp)*Sqr(Qq)*Sqr(Qs) + 1620*Quad(gp)*
      Sqr(Qd)*Sqr(Qu) + 1620*Quad(gp)*Sqr(QHd)*Sqr(Qu) + 1620*Quad(gp)*Sqr(Qq)*
      Sqr(Qu) + 540*Quad(gp)*Sqr(Qd)*Sqr(Qv) + 540*Quad(gp)*Sqr(QHd)*Sqr(Qv) +
      540*Quad(gp)*Sqr(Qq)*Sqr(Qv) - 270*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + 0.2
      *(-45*traceYdAdjYd - 15*traceYeAdjYe - 15*AbsSqr(Lambdax) + 4*Sqr(g1) +
      30*Sqr(g2) - 10*Sqr(gp)*Sqr(Qd) + 30*Sqr(gp)*Sqr(QHd) + 10*Sqr(gp)*Sqr(Qq
      ))*(Yd*Yd.adjoint()*Yd) + 0.2*(-15*traceYuAdjYu - 5*traceYvAdjYv - 5*
      AbsSqr(Lambdax) + 4*Sqr(g1) + 10*Sqr(gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Qq) +
      10*Sqr(gp)*Sqr(Qu))*(Yd*Yu.adjoint()*Yu) - 4*(Yd*Yd.adjoint()*Yd*Yd.
      adjoint()*Yd) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) - 2*(Yd*Yu.adjoint
      ()*Yu*Yu.adjoint()*Yu))).real();


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

/**
 * Calculates the 4-loop beta function of Yd.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yd_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

/**
 * Calculates the 5-loop beta function of Yd.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yd_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = ZEROMATRIX(3,3);


   return beta_Yd;
}

} // namespace flexiblesusy
