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

   beta_Yu = (-0.06666666666666667*Yu*(-45*traceYuAdjYu - 15*traceYvAdjYv - 15*
      AbsSqr(Lambdax) + 13*Sqr(g1) + 45*Sqr(g2) + 80*Sqr(g3) + 30*Sqr(gp)*Sqr(
      QHu) + 30*Sqr(gp)*Sqr(Qq) + 30*Sqr(gp)*Sqr(Qu)) + Yu*Yd.adjoint()*Yd + 3*
      (Yu*Yu.adjoint()*Yu)).real();


   return oneLoop * beta_Yu;
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
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;


   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (0.0022222222222222222*Yu*(-1350*traceYdAdjYuYuAdjYd - 4050*
      traceYuAdjYuYuAdjYu - 450*traceYvAdjYvTpYeconjYe - 1350*
      traceYvAdjYvYvAdjYv - 1350*traceYdAdjYd*AbsSqr(Lambdax) - 450*
      traceYeAdjYe*AbsSqr(Lambdax) + 2743*Quad(g1) + 3375*Quad(g2) - 800*Quad(
      g3) + 3600*Quad(gp)*Quad(QHu) + 18000*Quad(gp)*Quad(Qq) + 9900*Quad(gp)*
      Quad(Qu) + 360*traceYuAdjYu*Sqr(g1) + 450*Sqr(g1)*Sqr(g2) + 7200*
      traceYuAdjYu*Sqr(g3) + 1360*Sqr(g1)*Sqr(g3) + 3600*Sqr(g2)*Sqr(g3) + 1620
      *Qd*QHu*Sqr(g1)*Sqr(gp) + 1620*Qe*QHu*Sqr(g1)*Sqr(gp) - 540*QHd*QHu*Sqr(
      g1)*Sqr(gp) - 1620*QHu*Ql*Sqr(g1)*Sqr(gp) + 540*Qd*Qq*Sqr(g1)*Sqr(gp) +
      540*Qe*Qq*Sqr(g1)*Sqr(gp) - 180*QHd*Qq*Sqr(g1)*Sqr(gp) + 1800*QHu*Qq*Sqr(
      g1)*Sqr(gp) - 540*Ql*Qq*Sqr(g1)*Sqr(gp) - 2160*Qd*Qu*Sqr(g1)*Sqr(gp) -
      2160*Qe*Qu*Sqr(g1)*Sqr(gp) + 720*QHd*Qu*Sqr(g1)*Sqr(gp) - 3960*QHu*Qu*Sqr
      (g1)*Sqr(gp) + 2160*Ql*Qu*Sqr(g1)*Sqr(gp) - 3240*Qq*Qu*Sqr(g1)*Sqr(gp) +
      900*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) - 2700*traceYuAdjYu*Sqr(gp)*Sqr(QHu)
      - 900*traceYvAdjYv*Sqr(gp)*Sqr(QHu) - 900*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu
      ) + 1080*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 2700*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 8100*
      Quad(gp)*Sqr(Qd)*Sqr(QHu) + 2700*Quad(gp)*Sqr(Qe)*Sqr(QHu) + 1800*Quad(gp
      )*Sqr(QHd)*Sqr(QHu) + 900*traceYvAdjYv*Sqr(gp)*Sqr(Ql) + 5400*Quad(gp)*
      Sqr(QHu)*Sqr(Ql) + 2700*traceYuAdjYu*Sqr(gp)*Sqr(Qq) + 600*Sqr(g1)*Sqr(gp
      )*Sqr(Qq) + 2700*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 4800*Sqr(g3)*Sqr(gp)*Sqr(Qq) +
      8100*Quad(gp)*Sqr(Qd)*Sqr(Qq) + 2700*Quad(gp)*Sqr(Qe)*Sqr(Qq) + 1800*Quad
      (gp)*Sqr(QHd)*Sqr(Qq) + 18000*Quad(gp)*Sqr(QHu)*Sqr(Qq) + 5400*Quad(gp)*
      Sqr(Ql)*Sqr(Qq) + 900*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 900*Quad(gp)*Sqr(
      QHu)*Sqr(Qs) + 900*Quad(gp)*Sqr(Qq)*Sqr(Qs) + 2700*traceYuAdjYu*Sqr(gp)*
      Sqr(Qu) + 5280*Sqr(g1)*Sqr(gp)*Sqr(Qu) + 4800*Sqr(g3)*Sqr(gp)*Sqr(Qu) +
      8100*Quad(gp)*Sqr(Qd)*Sqr(Qu) + 2700*Quad(gp)*Sqr(Qe)*Sqr(Qu) + 1800*Quad
      (gp)*Sqr(QHd)*Sqr(Qu) + 9900*Quad(gp)*Sqr(QHu)*Sqr(Qu) + 5400*Quad(gp)*
      Sqr(Ql)*Sqr(Qu) + 24300*Quad(gp)*Sqr(Qq)*Sqr(Qu) + 900*Quad(gp)*Sqr(Qs)*
      Sqr(Qu) + 900*traceYvAdjYv*Sqr(gp)*Sqr(Qv) + 2700*Quad(gp)*Sqr(QHu)*Sqr(
      Qv) + 2700*Quad(gp)*Sqr(Qq)*Sqr(Qv) + 2700*Quad(gp)*Sqr(Qu)*Sqr(Qv) -
      1350*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + 0.2*(-15*traceYdAdjYd - 5*
      traceYeAdjYe - 5*AbsSqr(Lambdax) + 2*Sqr(g1) + 10*Sqr(gp)*Sqr(Qd) + 10*
      Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Qq))*(Yu*Yd.adjoint()*Yd) + 0.2*(-45*
      traceYuAdjYu - 15*traceYvAdjYv - 15*AbsSqr(Lambdax) + 2*Sqr(g1) + 30*Sqr(
      g2) + 30*Sqr(gp)*Sqr(QHu) + 10*Sqr(gp)*Sqr(Qq) - 10*Sqr(gp)*Sqr(Qu))*(Yu*
      Yu.adjoint()*Yu) - 2*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(Yu*Yd.
      adjoint()*Yd*Yu.adjoint()*Yu) - 4*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*Yu)).
      real();


   return twoLoop * beta_Yu;
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


   return threeLoop * beta_Yu;
}

/**
 * Calculates the 4-loop beta function of Yu.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yu_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return fourLoop * beta_Yu;
}

/**
 * Calculates the 5-loop beta function of Yu.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yu_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = ZEROMATRIX(3,3);


   return fiveLoop * beta_Yu;
}

} // namespace flexiblesusy
