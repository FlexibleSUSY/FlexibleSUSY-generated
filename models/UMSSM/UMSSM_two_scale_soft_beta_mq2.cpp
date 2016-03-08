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

// File generated at Tue 8 Mar 2016 18:00:32

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
      0.2581988897471611*g1*Tr11*UNITMATRIX(3) + 2*gp*Qq*Tr14*UNITMATRIX(3) -
      0.13333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) - 6*AbsSqr(MassWB
      )*Sqr(g2)*UNITMATRIX(3) - 10.666666666666666*AbsSqr(MassG)*Sqr(g3)*
      UNITMATRIX(3) - 8*AbsSqr(MassU)*Sqr(gp)*Sqr(Qq)*UNITMATRIX(3))).real();


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

   const Eigen::Matrix<double,3,3> beta_mq2_1 = (0.2*twoLoop*(2*(-15*
      traceconjTYdTpTYd - 5*traceconjTYeTpTYe - 15*tracemd2YdAdjYd - 5*
      traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*tracemq2AdjYdYd - 30*mHd2*
      traceYdAdjYd - 10*mHd2*traceYeAdjYe - 10*mHd2*AbsSqr(Lambdax) - 5*mHu2*
      AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) + 2*mHd2*Sqr(g1) + 4*AbsSqr(MassB
      )*Sqr(g1) + 10*mHd2*Sqr(gp)*Sqr(Qd) + 10*mHd2*Sqr(gp)*Sqr(QHd) + 20*
      AbsSqr(MassU)*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq)) - 10*mHd2*Sqr(gp)*
      Sqr(Qq))*(Yd.adjoint()*Yd) - 2*(2*Conj(MassB)*Sqr(g1) + 5*(3*
      traceconjTYdTpYd + traceconjTYeTpYe + Conj(TLambdax)*Lambdax + 2*Conj(
      MassU)*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq))))*(Yd.adjoint()*TYd) - 30*
      traceconjTYuTpTYu*(Yu.adjoint()*Yu) - 10*traceconjTYvTpTYv*(Yu.adjoint()*
      Yu) - 30*tracemq2AdjYuYu*(Yu.adjoint()*Yu) - 30*tracemu2YuAdjYu*(
      Yu.adjoint()*Yu) - 60*mHu2*traceYuAdjYu*(Yu.adjoint()*Yu) - 20*mHu2*
      traceYvAdjYv*(Yu.adjoint()*Yu) - 10*traceYvAdjYvconjml2*(Yu.adjoint()*Yu)
      - 10*traceYvconjmvR2AdjYv*(Yu.adjoint()*Yu) - 10*mHd2*AbsSqr(Lambdax)*(
      Yu.adjoint()*Yu) - 20*mHu2*AbsSqr(Lambdax)*(Yu.adjoint()*Yu) - 10*ms2*
      AbsSqr(Lambdax)*(Yu.adjoint()*Yu) + 8*mHu2*Sqr(g1)*(Yu.adjoint()*Yu) + 16
      *AbsSqr(MassB)*Sqr(g1)*(Yu.adjoint()*Yu) + 20*mHu2*Sqr(gp)*Sqr(QHu)*(
      Yu.adjoint()*Yu) + 40*AbsSqr(MassU)*Sqr(gp)*Sqr(QHu)*(Yu.adjoint()*Yu) -
      20*mHu2*Sqr(gp)*Sqr(Qq)*(Yu.adjoint()*Yu) - 40*AbsSqr(MassU)*Sqr(gp)*Sqr(
      Qq)*(Yu.adjoint()*Yu) + 20*mHu2*Sqr(gp)*Sqr(Qu)*(Yu.adjoint()*Yu) + 40*
      AbsSqr(MassU)*Sqr(gp)*Sqr(Qu)*(Yu.adjoint()*Yu) - 30*traceconjTYuTpYu*(
      Yu.adjoint()*TYu) - 10*traceconjTYvTpYv*(Yu.adjoint()*TYu) - 10*Conj(
      TLambdax)*Lambdax*(Yu.adjoint()*TYu) - 8*Conj(MassB)*Sqr(g1)*(Yu.adjoint(
      )*TYu) - 20*Conj(MassU)*Sqr(gp)*Sqr(QHu)*(Yu.adjoint()*TYu) + 20*Conj(
      MassU)*Sqr(gp)*Sqr(Qq)*(Yu.adjoint()*TYu) - 20*Conj(MassU)*Sqr(gp)*Sqr(Qu
      )*(Yu.adjoint()*TYu) - 30*traceAdjYdTYd*((TYd).adjoint()*Yd) - 10*
      traceAdjYeTYe*((TYd).adjoint()*Yd) - 4*MassB*Sqr(g1)*((TYd).adjoint()*Yd)
      - 20*MassU*Sqr(gp)*Sqr(Qd)*((TYd).adjoint()*Yd) - 20*MassU*Sqr(gp)*Sqr(
      QHd)*((TYd).adjoint()*Yd) + 20*MassU*Sqr(gp)*Sqr(Qq)*((TYd).adjoint()*Yd)
      - 30*traceYdAdjYd*((TYd).adjoint()*TYd) - 10*traceYeAdjYe*((TYd).adjoint
      ()*TYd) - 10*AbsSqr(Lambdax)*((TYd).adjoint()*TYd) + 4*Sqr(g1)*((TYd)
      .adjoint()*TYd) + 20*Sqr(gp)*Sqr(Qd)*((TYd).adjoint()*TYd) + 20*Sqr(gp)*
      Sqr(QHd)*((TYd).adjoint()*TYd) - 20*Sqr(gp)*Sqr(Qq)*((TYd).adjoint()*TYd)
      - 30*traceAdjYuTYu*((TYu).adjoint()*Yu) - 10*traceAdjYvTYv*((TYu)
      .adjoint()*Yu) - 8*MassB*Sqr(g1)*((TYu).adjoint()*Yu) - 20*MassU*Sqr(gp)*
      Sqr(QHu)*((TYu).adjoint()*Yu) + 20*MassU*Sqr(gp)*Sqr(Qq)*((TYu).adjoint()
      *Yu) - 20*MassU*Sqr(gp)*Sqr(Qu)*((TYu).adjoint()*Yu) - 30*traceYuAdjYu*((
      TYu).adjoint()*TYu) - 10*traceYvAdjYv*((TYu).adjoint()*TYu) - 10*AbsSqr(
      Lambdax)*((TYu).adjoint()*TYu) + 8*Sqr(g1)*((TYu).adjoint()*TYu) + 20*Sqr
      (gp)*Sqr(QHu)*((TYu).adjoint()*TYu) - 20*Sqr(gp)*Sqr(Qq)*((TYu).adjoint()
      *TYu) + 20*Sqr(gp)*Sqr(Qu)*((TYu).adjoint()*TYu) - 15*traceYdAdjYd*(mq2*
      Yd.adjoint()*Yd) - 5*traceYeAdjYe*(mq2*Yd.adjoint()*Yd) - 5*AbsSqr(
      Lambdax)*(mq2*Yd.adjoint()*Yd) + 2*Sqr(g1)*(mq2*Yd.adjoint()*Yd) + 10*Sqr
      (gp)*Sqr(Qd)*(mq2*Yd.adjoint()*Yd) + 10*Sqr(gp)*Sqr(QHd)*(mq2*Yd.adjoint(
      )*Yd) - 10*Sqr(gp)*Sqr(Qq)*(mq2*Yd.adjoint()*Yd) - 15*traceYuAdjYu*(mq2*
      Yu.adjoint()*Yu) - 5*traceYvAdjYv*(mq2*Yu.adjoint()*Yu) - 5*AbsSqr(
      Lambdax)*(mq2*Yu.adjoint()*Yu) + 4*Sqr(g1)*(mq2*Yu.adjoint()*Yu) + 10*Sqr
      (gp)*Sqr(QHu)*(mq2*Yu.adjoint()*Yu) - 10*Sqr(gp)*Sqr(Qq)*(mq2*Yu.adjoint(
      )*Yu) + 10*Sqr(gp)*Sqr(Qu)*(mq2*Yu.adjoint()*Yu) - 30*traceYdAdjYd*(
      Yd.adjoint()*md2*Yd) - 10*traceYeAdjYe*(Yd.adjoint()*md2*Yd) - 10*AbsSqr(
      Lambdax)*(Yd.adjoint()*md2*Yd) + 4*Sqr(g1)*(Yd.adjoint()*md2*Yd) + 20*Sqr
      (gp)*Sqr(Qd)*(Yd.adjoint()*md2*Yd) + 20*Sqr(gp)*Sqr(QHd)*(Yd.adjoint()*
      md2*Yd) - 20*Sqr(gp)*Sqr(Qq)*(Yd.adjoint()*md2*Yd) + 2*Sqr(g1)*(
      Yd.adjoint()*Yd*mq2))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_mq2_2 = (twoLoop*UNITMATRIX(3)*(
      -2*AbsSqr(TLambdax)*(Yd.adjoint()*Yd) - 2*AbsSqr(TLambdax)*(Yu.adjoint()*
      Yu) - 2*Conj(Lambdax)*TLambdax*((TYd).adjoint()*Yd) - 2*Conj(Lambdax)*
      TLambdax*((TYu).adjoint()*Yu) + (-3*traceYdAdjYd - traceYeAdjYe - AbsSqr(
      Lambdax) + 2*Sqr(gp)*(Sqr(Qd) + Sqr(QHd) - Sqr(Qq)))*(Yd.adjoint()*Yd*mq2
      ) + 0.4*(-5*AbsSqr(Lambdax) + 4*Sqr(g1) + 5*(-3*traceYuAdjYu -
      traceYvAdjYv + 2*Sqr(gp)*(Sqr(QHu) - Sqr(Qq) + Sqr(Qu))))*(Yu.adjoint()*
      mu2*Yu) - 3*traceYuAdjYu*(Yu.adjoint()*Yu*mq2) - traceYvAdjYv*(Yu.adjoint
      ()*Yu*mq2) - AbsSqr(Lambdax)*(Yu.adjoint()*Yu*mq2) + 0.8*Sqr(g1)*(
      Yu.adjoint()*Yu*mq2) + 2*Sqr(gp)*Sqr(QHu)*(Yu.adjoint()*Yu*mq2) - 2*Sqr(
      gp)*Sqr(Qq)*(Yu.adjoint()*Yu*mq2) + 2*Sqr(gp)*Sqr(Qu)*(Yu.adjoint()*Yu*
      mq2) - 8*mHd2*(Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 4*(Yd.adjoint()*Yd*(TYd
      ).adjoint()*TYd) - 4*(Yd.adjoint()*TYd*(TYd).adjoint()*Yd) - 8*mHu2*(
      Yu.adjoint()*Yu*Yu.adjoint()*Yu) - 4*(Yu.adjoint()*Yu*(TYu).adjoint()*TYu
      ) - 4*(Yu.adjoint()*TYu*(TYu).adjoint()*Yu) - 4*((TYd).adjoint()*Yd*
      Yd.adjoint()*TYd) - 4*((TYd).adjoint()*TYd*Yd.adjoint()*Yd) - 4*((TYu)
      .adjoint()*Yu*Yu.adjoint()*TYu) - 4*((TYu).adjoint()*TYu*Yu.adjoint()*Yu)
      - 2*(mq2*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(mq2*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu) - 4*(Yd.adjoint()*md2*Yd*Yd.adjoint()*Yd) - 4*(
      Yd.adjoint()*Yd*mq2*Yd.adjoint()*Yd) - 4*(Yd.adjoint()*Yd*Yd.adjoint()*
      md2*Yd) - 2*(Yd.adjoint()*Yd*Yd.adjoint()*Yd*mq2) - 4*(Yu.adjoint()*mu2*
      Yu*Yu.adjoint()*Yu) - 4*(Yu.adjoint()*Yu*mq2*Yu.adjoint()*Yu) - 4*(
      Yu.adjoint()*Yu*Yu.adjoint()*mu2*Yu) - 2*(Yu.adjoint()*Yu*Yu.adjoint()*Yu
      *mq2) + 6*Power(g2,4)*Tr22*UNITMATRIX(3) + 10.666666666666666*Power(g3,4)
      *Tr23*UNITMATRIX(3) + 1.0327955589886444*g1*gp*Qq*Tr2U114*UNITMATRIX(3) +
      1.0327955589886444*g1*gp*Qq*Tr2U141*UNITMATRIX(3) + 1.0327955589886444*
      g1*Tr31*UNITMATRIX(3) + 8*gp*Qq*Tr34*UNITMATRIX(3) + 2.6533333333333333*
      Power(g1,4)*AbsSqr(MassB)*UNITMATRIX(3) - 42.666666666666664*Power(g3,4)*
      AbsSqr(MassG)*UNITMATRIX(3) + 480*Power(gp,4)*Power(Qq,4)*AbsSqr(MassU)*
      UNITMATRIX(3) + 0.13333333333333333*Tr2U111*Sqr(g1)*UNITMATRIX(3) + 0.4*
      AbsSqr(MassB)*Sqr(g1)*Sqr(g2)*UNITMATRIX(3) + 0.2*MassWB*Conj(MassB)*Sqr(
      g1)*Sqr(g2)*UNITMATRIX(3) + 0.7111111111111111*AbsSqr(MassB)*Sqr(g1)*Sqr(
      g3)*UNITMATRIX(3) + 0.7111111111111111*AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 0.35555555555555557*MassG*Conj(MassB)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 0.35555555555555557*MassB*Conj(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 32*AbsSqr(MassG)*Sqr(g2)*Sqr(g3)*UNITMATRIX(3) + 16*
      MassWB*Conj(MassG)*Sqr(g2)*Sqr(g3)*UNITMATRIX(3) + 4.8*Qd*Qq*AbsSqr(MassB
      )*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 4.8*Qe*Qq*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)
      *UNITMATRIX(3) - 1.6*QHd*Qq*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) +
      1.6*QHu*Qq*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 4.8*Ql*Qq*
      AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 9.6*Qq*Qu*AbsSqr(MassB)*Sqr
      (g1)*Sqr(gp)*UNITMATRIX(3) + 4.8*Qd*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*
      UNITMATRIX(3) + 4.8*Qe*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) -
      1.6*QHd*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 1.6*QHu*Qq*
      AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 4.8*Ql*Qq*AbsSqr(MassU)*Sqr
      (g1)*Sqr(gp)*UNITMATRIX(3) + 2.4*MassU*Qd*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp)*
      UNITMATRIX(3) + 2.4*MassU*Qe*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3)
      - 0.8*MassU*QHd*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 0.8*MassU
      *QHu*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 2.4*MassU*Ql*Qq*Conj(
      MassB)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 4.8*MassU*Qq*Qu*Conj(MassB)*Sqr(g1
      )*Sqr(gp)*UNITMATRIX(3) + 2.4*MassB*Qd*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp)*
      UNITMATRIX(3) + 2.4*MassB*Qe*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3)
      - 0.8*MassB*QHd*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 0.8*MassB
      *QHu*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) - 2.4*MassB*Ql*Qq*Conj(
      MassU)*Sqr(g1)*Sqr(gp)*UNITMATRIX(3) + 8*Tr2U144*Sqr(gp)*Sqr(Qq)*
      UNITMATRIX(3) + 5.333333333333333*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*Sqr(Qq)*
      UNITMATRIX(3) + 5.333333333333333*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*Sqr(Qq)*
      UNITMATRIX(3) + 2.6666666666666665*MassU*Conj(MassB)*Sqr(g1)*Sqr(gp)*Sqr(
      Qq)*UNITMATRIX(3) + 2.6666666666666665*MassB*Conj(MassU)*Sqr(g1)*Sqr(gp)*
      Sqr(Qq)*UNITMATRIX(3) + 24*AbsSqr(MassU)*Sqr(g2)*Sqr(gp)*Sqr(Qq)*
      UNITMATRIX(3) + 12*MassWB*Conj(MassU)*Sqr(g2)*Sqr(gp)*Sqr(Qq)*UNITMATRIX(
      3) + 42.666666666666664*AbsSqr(MassG)*Sqr(g3)*Sqr(gp)*Sqr(Qq)*UNITMATRIX(
      3) + 42.666666666666664*AbsSqr(MassU)*Sqr(g3)*Sqr(gp)*Sqr(Qq)*UNITMATRIX(
      3) + 21.333333333333332*MassU*Conj(MassG)*Sqr(g3)*Sqr(gp)*Sqr(Qq)*
      UNITMATRIX(3) + 21.333333333333332*MassG*Conj(MassU)*Sqr(g3)*Sqr(gp)*Sqr(
      Qq)*UNITMATRIX(3) + 216*Power(gp,4)*AbsSqr(MassU)*Sqr(Qd)*Sqr(Qq)*
      UNITMATRIX(3) + 72*Power(gp,4)*AbsSqr(MassU)*Sqr(Qe)*Sqr(Qq)*UNITMATRIX(3
      ) + 48*Power(gp,4)*AbsSqr(MassU)*Sqr(QHd)*Sqr(Qq)*UNITMATRIX(3) + 48*
      Power(gp,4)*AbsSqr(MassU)*Sqr(QHu)*Sqr(Qq)*UNITMATRIX(3) + 144*Power(gp,4
      )*AbsSqr(MassU)*Sqr(Ql)*Sqr(Qq)*UNITMATRIX(3))).real();
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
