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

// File generated at Wed 12 Apr 2017 12:26:52

#include "UMSSM_two_scale_soft_parameters.hpp"
#include "wrappers.hpp"

namespace flexiblesusy {

#define INPUT(parameter) input.parameter
#define TRACE_STRUCT soft_traces

namespace {

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator+(double n, const Eigen::MatrixBase<Derived>& m)
{
   return m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(const Eigen::MatrixBase<Derived>& m, double n)
{
   return m - Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

template <typename Derived>
typename Eigen::MatrixBase<Derived>::PlainObject operator-(double n, const Eigen::MatrixBase<Derived>& m)
{
   return - m + Eigen::Matrix<double,
                            Eigen::MatrixBase<Derived>::RowsAtCompileTime,
                            Eigen::MatrixBase<Derived>::ColsAtCompileTime>::Identity() * n;
}

} // anonymous namespace

/**
 * Calculates the one-loop beta function of mq2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mq2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qq = INPUT(Qq);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = (oneOver16PiSqr*(2*mHd2*(Yd.adjoint()*Yd) + 2*mHu2*(
      Yu.adjoint()*Yu) + 2*((TYd).adjoint()*TYd) + 2*((TYu).adjoint()*TYu) +
      mq2*Yd.adjoint()*Yd + mq2*Yu.adjoint()*Yu + 2*(Yd.adjoint()*md2*Yd) +
      Yd.adjoint()*Yd*mq2 + 2*(Yu.adjoint()*mu2*Yu) + Yu.adjoint()*Yu*mq2 +
      0.06666666666666667*(3.872983346207417*g1*Tr11 + 30*gp*Qq*Tr14 - 2*AbsSqr
      (MassB)*Sqr(g1) - 90*AbsSqr(MassWB)*Sqr(g2) - 160*AbsSqr(MassG)*Sqr(g3) -
      120*AbsSqr(MassU)*Sqr(gp)*Sqr(Qq))*UNITMATRIX(3))).real();


   return beta_mq2;
}

/**
 * Calculates the two-loop beta function of mq2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mq2_two_loop(const Soft_traces& soft_traces) const
{
   const auto Qq = INPUT(Qq);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mq2;

   const Eigen::Matrix<double,3,3> beta_mq2_1 = ((0.4*twoLoop*(-15*
      traceconjTYdTpTYd - 5*traceconjTYeTpTYe - 15*tracemd2YdAdjYd - 5*
      traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*tracemq2AdjYdYd - 30*mHd2*
      traceYdAdjYd - 10*mHd2*traceYeAdjYe - 10*mHd2*AbsSqr(Lambdax) - 5*mHu2*
      AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) + 2*mHd2*Sqr(g1) + 4*AbsSqr(MassB
      )*Sqr(g1) + 10*mHd2*Sqr(gp)*Sqr(Qd) + 10*mHd2*Sqr(gp)*Sqr(QHd) + 20*
      AbsSqr(MassU)*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq)) - 10*mHd2*Sqr(gp)*
      Sqr(Qq))*(Yd.adjoint()*Yd) - 0.4*twoLoop*(2*Conj(MassB)*Sqr(g1) + 5*(3*
      traceconjTYdTpYd + traceconjTYeTpYe + Conj(TLambdax)*Lambdax + 2*Conj(
      MassU)*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq))))*(Yd.adjoint()*TYd) + 0.4*
      twoLoop*(-15*traceconjTYuTpTYu - 5*traceconjTYvTpTYv - 15*tracemq2AdjYuYu
      - 15*tracemu2YuAdjYu - 30*mHu2*traceYuAdjYu - 10*mHu2*traceYvAdjYv - 5*
      traceYvAdjYvconjml2 - 5*traceYvconjmvR2AdjYv - 5*mHd2*AbsSqr(Lambdax) -
      10*mHu2*AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) + 4*mHu2*Sqr(g1) + 8*
      AbsSqr(MassB)*Sqr(g1) + 10*mHu2*Sqr(gp)*Sqr(QHu) - 10*mHu2*Sqr(gp)*Sqr(Qq
      ) + 10*mHu2*Sqr(gp)*Sqr(Qu) + 20*AbsSqr(MassU)*Sqr(gp)*(Sqr(QHu) - Sqr(Qq
      ) + Sqr(Qu)))*(Yu.adjoint()*Yu) - 0.4*twoLoop*(4*Conj(MassB)*Sqr(g1) + 5*
      (3*traceconjTYuTpYu + traceconjTYvTpYv + Conj(TLambdax)*Lambdax + 2*Conj(
      MassU)*Sqr(gp)*(Sqr(QHu) - Sqr(Qq) + Sqr(Qu))))*(Yu.adjoint()*TYu) - 0.4*
      twoLoop*(2*MassB*Sqr(g1) + 5*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassU*
      Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq))))*((TYd).adjoint()*Yd) + 0.4*
      twoLoop*(-5*AbsSqr(Lambdax) + 2*Sqr(g1) + 5*(-3*traceYdAdjYd -
      traceYeAdjYe + 2*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq))))*((TYd).adjoint(
      )*TYd) - 0.4*twoLoop*(4*MassB*Sqr(g1) + 5*(3*traceAdjYuTYu +
      traceAdjYvTYv + 2*MassU*Sqr(gp)*(Sqr(QHu) - Sqr(Qq) + Sqr(Qu))))*((TYu)
      .adjoint()*Yu) + 0.4*twoLoop*(-5*AbsSqr(Lambdax) + 4*Sqr(g1) + 5*(-3*
      traceYuAdjYu - traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHu) - Sqr(Qq) + Sqr(Qu))))*
      ((TYu).adjoint()*TYu) + 0.2*twoLoop*(-5*AbsSqr(Lambdax) + 2*Sqr(g1) + 5*(
      -3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq))
      ))*(mq2*Yd.adjoint()*Yd) + 0.2*twoLoop*(-5*AbsSqr(Lambdax) + 4*Sqr(g1) +
      5*(-3*traceYuAdjYu - traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHu) - Sqr(Qq) + Sqr(
      Qu))))*(mq2*Yu.adjoint()*Yu) + 0.4*twoLoop*(-5*AbsSqr(Lambdax) + 2*Sqr(g1
      ) + 5*(-3*traceYdAdjYd - traceYeAdjYe + 2*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) -
      Sqr(Qq))))*(Yd.adjoint()*md2*Yd) + 0.4*twoLoop*Sqr(g1)*(Yd.adjoint()*Yd*
      mq2))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_mq2_2 = (UNITMATRIX(3)*(-2*
      twoLoop*AbsSqr(TLambdax)*(Yd.adjoint()*Yd) - 2*twoLoop*AbsSqr(TLambdax)*(
      Yu.adjoint()*Yu) - 2*twoLoop*Conj(Lambdax)*TLambdax*((TYd).adjoint()*Yd)
      - 2*twoLoop*Conj(Lambdax)*TLambdax*((TYu).adjoint()*Yu) + twoLoop*(-3*
      traceYdAdjYd - traceYeAdjYe - AbsSqr(Lambdax) + 2*Sqr(gp)*(Sqr(Qd) + Sqr(
      QHd) - Sqr(Qq)))*(Yd.adjoint()*Yd*mq2) + 0.4*twoLoop*(-5*AbsSqr(Lambdax)
      + 4*Sqr(g1) + 5*(-3*traceYuAdjYu - traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHu) -
      Sqr(Qq) + Sqr(Qu))))*(Yu.adjoint()*mu2*Yu) + 0.2*twoLoop*(-5*AbsSqr(
      Lambdax) + 4*Sqr(g1) + 5*(-3*traceYuAdjYu - traceYvAdjYv + 2*Sqr(gp)*(Sqr
      (QHu) - Sqr(Qq) + Sqr(Qu))))*(Yu.adjoint()*Yu*mq2) - 8*mHd2*twoLoop*(
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*twoLoop*(Yd.adjoint()*Yd*(TYd)
      .adjoint()*TYd) - 4*twoLoop*(Yd.adjoint()*TYd*(TYd).adjoint()*Yd) - 8*
      mHu2*twoLoop*(Yu.adjoint()*Yu*Yu.adjoint()*Yu) - 4*twoLoop*(Yu.adjoint()*
      Yu*(TYu).adjoint()*TYu) - 4*twoLoop*(Yu.adjoint()*TYu*(TYu).adjoint()*Yu)
      - 4*twoLoop*((TYd).adjoint()*Yd*Yd.adjoint()*TYd) - 4*twoLoop*((TYd)
      .adjoint()*TYd*Yd.adjoint()*Yd) - 4*twoLoop*((TYu).adjoint()*Yu*
      Yu.adjoint()*TYu) - 4*twoLoop*((TYu).adjoint()*TYu*Yu.adjoint()*Yu) - 2*
      twoLoop*(mq2*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*twoLoop*(mq2*Yu.adjoint
      ()*Yu*Yu.adjoint()*Yu) - 4*twoLoop*(Yd.adjoint()*md2*Yd*Yd.adjoint()*Yd)
      - 4*twoLoop*(Yd.adjoint()*Yd*mq2*Yd.adjoint()*Yd) - 4*twoLoop*(Yd.adjoint
      ()*Yd*Yd.adjoint()*md2*Yd) - 2*twoLoop*(Yd.adjoint()*Yd*Yd.adjoint()*Yd*
      mq2) - 4*twoLoop*(Yu.adjoint()*mu2*Yu*Yu.adjoint()*Yu) - 4*twoLoop*(
      Yu.adjoint()*Yu*mq2*Yu.adjoint()*Yu) - 4*twoLoop*(Yu.adjoint()*Yu*
      Yu.adjoint()*mu2*Yu) - 2*twoLoop*(Yu.adjoint()*Yu*Yu.adjoint()*Yu*mq2) +
      0.0044444444444444444*twoLoop*(Conj(MassB)*Sqr(g1)*(597*MassB*Sqr(g1) + 5
      *(9*(2*MassB + MassWB)*Sqr(g2) + 4*(4*(2*MassB + MassG)*Sqr(g3) + 3*(2*
      MassB + MassU)*Qq*(9*Qd + 9*Qe - 3*QHd + 3*QHu - 9*Ql + 10*Qq - 18*Qu)*
      Sqr(gp)))) + 10*(8*Conj(MassG)*Sqr(g3)*((MassB + 2*MassG)*Sqr(g1) + 15*(3
      *(2*MassG + MassWB)*Sqr(g2) - 8*MassG*Sqr(g3) + 4*(2*MassG + MassU)*Sqr(
      gp)*Sqr(Qq))) + 3*(45*Power(g2,4)*Tr22 + 80*Power(g3,4)*Tr23 +
      7.745966692414834*g1*gp*Qq*Tr2U114 + 7.745966692414834*g1*gp*Qq*Tr2U141 +
      7.745966692414834*g1*Tr31 + 60*gp*Qq*Tr34 + Tr2U111*Sqr(g1) + 60*Tr2U144
      *Sqr(gp)*Sqr(Qq) + 2*Qq*Conj(MassU)*Sqr(gp)*((MassB + 2*MassU)*(9*Qd + 9*
      Qe - 3*QHd + 3*QHu - 9*Ql + 10*Qq)*Sqr(g1) + 5*Qq*(16*(MassG + 2*MassU)*
      Sqr(g3) + 9*((2*MassU + MassWB)*Sqr(g2) + 2*MassU*Sqr(gp)*(9*Sqr(Qd) + 3*
      Sqr(Qe) + 2*(Sqr(QHd) + Sqr(QHu) + 3*Sqr(Ql) + 10*Sqr(Qq)))))))))*
      UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_mq2_3 = (0.2*twoLoop*(Conj(MassWB
      )*Sqr(g2)*((MassB + 2*MassWB)*Sqr(g1) + 5*(16*(MassG + 2*MassWB)*Sqr(g3)
      + 3*(11*MassWB*Sqr(g2) + 4*(MassU + 2*MassWB)*Sqr(gp)*Sqr(Qq)))) + 24*Qq*
      Conj(MassU)*Sqr(gp)*(-((MassB + 2*MassU)*Qu*Sqr(g1)) + 5*MassU*Qq*Sqr(gp)
      *(Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv))))*UNITMATRIX(3)).real();

   beta_mq2 = beta_mq2_1 + beta_mq2_2 + beta_mq2_3;


   return beta_mq2;
}

/**
 * Calculates the three-loop beta function of mq2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mq2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = ZEROMATRIX(3,3);


   return beta_mq2;
}

} // namespace flexiblesusy
