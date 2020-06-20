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
 * Calculates the 1-loop beta function of Yv.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_1_loop(const Susy_traces& susy_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qv = INPUT(Qv);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;


   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = (-0.2*Yv*(-15*traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) +
      3*Sqr(g1) + 15*Sqr(g2) + 10*Sqr(gp)*Sqr(QHu) + 10*Sqr(gp)*Sqr(Ql) + 10*
      Sqr(gp)*Sqr(Qv)) + 3*(Yv*Yv.adjoint()*Yv) + Ye.transpose()*Ye.conjugate()
      *Yv).real();


   return oneLoop * beta_Yv;
}

/**
 * Calculates the 2-loop beta function of Yv.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_2_loop(const Susy_traces& susy_traces) const
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


   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = (0.02*Yv*(-150*traceYdAdjYuYuAdjYd - 450*traceYuAdjYuYuAdjYu - 50*
      traceYvAdjYvTpYeconjYe - 150*traceYvAdjYvYvAdjYv - 150*traceYdAdjYd*
      AbsSqr(Lambdax) - 50*traceYeAdjYe*AbsSqr(Lambdax) + 207*Quad(g1) + 375*
      Quad(g2) + 400*Quad(gp)*Quad(QHu) + 800*Quad(gp)*Quad(Ql) + 500*Quad(gp)*
      Quad(Qv) + 40*traceYuAdjYu*Sqr(g1) + 90*Sqr(g1)*Sqr(g2) + 800*
      traceYuAdjYu*Sqr(g3) + 180*Qd*QHu*Sqr(g1)*Sqr(gp) + 180*Qe*QHu*Sqr(g1)*
      Sqr(gp) - 60*QHd*QHu*Sqr(g1)*Sqr(gp) - 180*Qd*Ql*Sqr(g1)*Sqr(gp) - 180*Qe
      *Ql*Sqr(g1)*Sqr(gp) + 60*QHd*Ql*Sqr(g1)*Sqr(gp) - 240*QHu*Ql*Sqr(g1)*Sqr(
      gp) + 180*QHu*Qq*Sqr(g1)*Sqr(gp) - 180*Ql*Qq*Sqr(g1)*Sqr(gp) - 360*QHu*Qu
      *Sqr(g1)*Sqr(gp) + 360*Ql*Qu*Sqr(g1)*Sqr(gp) + 100*AbsSqr(Lambdax)*Sqr(gp
      )*Sqr(QHd) - 300*traceYuAdjYu*Sqr(gp)*Sqr(QHu) - 100*traceYvAdjYv*Sqr(gp)
      *Sqr(QHu) - 100*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 120*Sqr(g1)*Sqr(gp)*
      Sqr(QHu) + 300*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 900*Quad(gp)*Sqr(Qd)*Sqr(QHu) +
      300*Quad(gp)*Sqr(Qe)*Sqr(QHu) + 200*Quad(gp)*Sqr(QHd)*Sqr(QHu) + 100*
      traceYvAdjYv*Sqr(gp)*Sqr(Ql) + 240*Sqr(g1)*Sqr(gp)*Sqr(Ql) + 300*Sqr(g2)*
      Sqr(gp)*Sqr(Ql) + 900*Quad(gp)*Sqr(Qd)*Sqr(Ql) + 300*Quad(gp)*Sqr(Qe)*Sqr
      (Ql) + 200*Quad(gp)*Sqr(QHd)*Sqr(Ql) + 800*Quad(gp)*Sqr(QHu)*Sqr(Ql) +
      300*traceYuAdjYu*Sqr(gp)*Sqr(Qq) + 1800*Quad(gp)*Sqr(QHu)*Sqr(Qq) + 1800*
      Quad(gp)*Sqr(Ql)*Sqr(Qq) + 100*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 100*Quad
      (gp)*Sqr(QHu)*Sqr(Qs) + 100*Quad(gp)*Sqr(Ql)*Sqr(Qs) + 300*traceYuAdjYu*
      Sqr(gp)*Sqr(Qu) + 900*Quad(gp)*Sqr(QHu)*Sqr(Qu) + 900*Quad(gp)*Sqr(Ql)*
      Sqr(Qu) + 100*traceYvAdjYv*Sqr(gp)*Sqr(Qv) + 900*Quad(gp)*Sqr(Qd)*Sqr(Qv)
      + 300*Quad(gp)*Sqr(Qe)*Sqr(Qv) + 200*Quad(gp)*Sqr(QHd)*Sqr(Qv) + 500*Quad
      (gp)*Sqr(QHu)*Sqr(Qv) + 900*Quad(gp)*Sqr(Ql)*Sqr(Qv) + 1800*Quad(gp)*Sqr(
      Qq)*Sqr(Qv) + 100*Quad(gp)*Sqr(Qs)*Sqr(Qv) + 900*Quad(gp)*Sqr(Qu)*Sqr(Qv)
      - 150*Sqr(Conj(Lambdax))*Sqr(Lambdax)) + 0.2*(-45*traceYuAdjYu - 15*
      traceYvAdjYv - 15*AbsSqr(Lambdax) + 6*Sqr(g1) + 30*Sqr(g2) + 30*Sqr(gp)*
      Sqr(QHu) + 10*Sqr(gp)*Sqr(Ql) - 10*Sqr(gp)*Sqr(Qv))*(Yv*Yv.adjoint()*Yv)
      + 0.2*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) + 6*Sqr(g1)
      + 10*Sqr(gp)*Sqr(Qe) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Ql))*(Ye.
      transpose()*Ye.conjugate()*Yv) - 4*(Yv*Yv.adjoint()*Yv*Yv.adjoint()*Yv) -
      2*(Yv*Yv.adjoint()*Ye.transpose()*Ye.conjugate()*Yv) - 2*(Ye.transpose()*
      Ye.conjugate()*Ye.transpose()*Ye.conjugate()*Yv)).real();


   return twoLoop * beta_Yv;
}

/**
 * Calculates the 3-loop beta function of Yv.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_3_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = ZEROMATRIX(3,3);


   return threeLoop * beta_Yv;
}

/**
 * Calculates the 4-loop beta function of Yv.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_4_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = ZEROMATRIX(3,3);


   return fourLoop * beta_Yv;
}

/**
 * Calculates the 5-loop beta function of Yv.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_susy_parameters::calc_beta_Yv_5_loop(const Susy_traces& susy_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_Yv;

   beta_Yv = ZEROMATRIX(3,3);


   return fiveLoop * beta_Yv;
}

} // namespace flexiblesusy
