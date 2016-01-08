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

// File generated at Fri 8 Jan 2016 15:12:17

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
 * Calculates the one-loop beta function of mHd2.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHd2;

   beta_mHd2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 + 2*gp*QHd*
      Tr14 + 6*traceconjTYdTpTYd + 2*traceconjTYeTpTYe + 6*tracemd2YdAdjYd + 2*
      traceme2YeAdjYe + 2*traceml2AdjYeYe + 6*tracemq2AdjYdYd + 6*mHd2*
      traceYdAdjYd + 2*mHd2*traceYeAdjYe + 2*mHd2*AbsSqr(Lambdax) + 2*mHu2*
      AbsSqr(Lambdax) + 2*ms2*AbsSqr(Lambdax) + 2*AbsSqr(TLambdax) - 1.2*AbsSqr
      (MassB)*Sqr(g1) - 6*AbsSqr(MassWB)*Sqr(g2) - 8*AbsSqr(MassU)*Sqr(gp)*Sqr(
      QHd)));


   return beta_mHd2;
}

/**
 * Calculates the two-loop beta function of mHd2.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qs = INPUT(Qs);
   const auto Qd = INPUT(Qd);
   const auto Qq = INPUT(Qq);
   const auto Qe = INPUT(Qe);
   const auto Ql = INPUT(Ql);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd =
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd =
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe =
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe =
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYeconjTYvTpTYvAdjYe =
      TRACE_STRUCT.traceYeconjTYvTpTYvAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;
   const double traceYvAdjYvTpTYeconjTYe =
      TRACE_STRUCT.traceYvAdjYvTpTYeconjTYe;
   const double traceAdjYeTYeconjTYvTpYv =
      TRACE_STRUCT.traceAdjYeTYeconjTYvTpYv;
   const double traceAdjYvTpYeconjTYeTYv =
      TRACE_STRUCT.traceAdjYvTpYeconjTYeTYv;
   const double tracemd2YdAdjYdYdAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe =
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe =
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double traceYvAdjYvconjml2TpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvconjml2TpYeconjYe;
   const double traceYvAdjYvTpYeconjme2conjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjme2conjYe;
   const double traceYvAdjYvTpYeconjYeconjml2 =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYeconjml2;
   const double traceYvconjmvR2AdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvconjmvR2AdjYvTpYeconjYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHd2;

   const double beta_mHd2_1 = Re(0.04*twoLoop*(3*Conj(MassB)*Sqr(g1)*(261
      *MassB*Sqr(g1) + 5*(3*(2*MassB + MassWB)*Sqr(g2) - 4*QHd*(MassU*(3*Qd + 3
      *Qe - 2*QHd + QHu - 3*Ql + 3*Qq - 6*Qu) + 2*MassB*(3*Qd + 3*Qe - 2*QHd +
      QHu - 3*Ql + 3*Qq - 6*Qu + 3*Qv))*Sqr(gp))) + 10*(15*Power(g2,4)*Tr22 -
      7.745966692414834*g1*(gp*QHd*(Tr2U114 + Tr2U141) + Tr31) + (3*Tr2U111 - 2
      *(traceconjTYdTpTYd - MassB*traceconjTYdTpYd - 3*traceconjTYeTpTYe + 3*
      MassB*traceconjTYeTpYe + tracemd2YdAdjYd - 3*traceme2YeAdjYe - 3*
      traceml2AdjYeYe + tracemq2AdjYdYd + mHd2*traceYdAdjYd - 3*mHd2*
      traceYeAdjYe))*Sqr(g1) - 5*(-4*gp*QHd*Tr34 + traceAdjYeTYeconjTYvTpYv +
      traceAdjYvTpYeconjTYeTYv + 18*tracemd2YdAdjYdYdAdjYd + 3*
      tracemd2YdAdjYuYuAdjYd + 6*traceme2YeAdjYeYeAdjYe + 6*
      traceml2AdjYeYeAdjYeYe + 18*tracemq2AdjYdYdAdjYdYd + 3*
      tracemq2AdjYdYdAdjYuYu + 3*tracemq2AdjYuYuAdjYdYd + 3*
      tracemu2YuAdjYdYdAdjYu + 18*traceYdAdjTYdTYdAdjYd + 3*
      traceYdAdjTYuTYuAdjYd + 18*traceYdAdjYdTYdAdjTYd + 18*mHd2*
      traceYdAdjYdYdAdjYd + 3*traceYdAdjYuTYuAdjTYd + 3*mHd2*
      traceYdAdjYuYuAdjYd + 3*mHu2*traceYdAdjYuYuAdjYd + 6*
      traceYeAdjTYeTYeAdjYe + 6*traceYeAdjYeTYeAdjTYe + 6*mHd2*
      traceYeAdjYeYeAdjYe + traceYeconjTYvTpTYvAdjYe + 3*traceYuAdjTYdTYdAdjYu
      + 3*traceYuAdjYdTYdAdjTYu + traceYvAdjYvconjml2TpYeconjYe +
      traceYvAdjYvTpTYeconjTYe + traceYvAdjYvTpYeconjme2conjYe + mHd2*
      traceYvAdjYvTpYeconjYe + mHu2*traceYvAdjYvTpYeconjYe +
      traceYvAdjYvTpYeconjYeconjml2 + traceYvconjmvR2AdjYvTpYeconjYe - 16*
      traceconjTYdTpTYd*Sqr(g3) + 16*MassG*traceconjTYdTpYd*Sqr(g3) - 16*
      tracemd2YdAdjYd*Sqr(g3) - 16*tracemq2AdjYdYd*Sqr(g3) - 16*mHd2*
      traceYdAdjYd*Sqr(g3) - 2*Sqr(gp)*(3*(traceconjTYdTpTYd - MassU*
      traceconjTYdTpYd + tracemd2YdAdjYd + tracemq2AdjYdYd + mHd2*traceYdAdjYd)
      *Sqr(Qd) + traceconjTYeTpTYe*Sqr(Qe) - MassU*traceconjTYeTpYe*Sqr(Qe) +
      traceme2YeAdjYe*Sqr(Qe) + traceml2AdjYeYe*Sqr(Qe) + mHd2*traceYeAdjYe*Sqr
      (Qe) + (2*Tr2U144 - 3*traceconjTYdTpTYd + 3*MassU*traceconjTYdTpYd -
      traceconjTYeTpTYe + MassU*traceconjTYeTpYe - 3*tracemd2YdAdjYd -
      traceme2YeAdjYe - traceml2AdjYeYe - 3*tracemq2AdjYdYd - 3*mHd2*
      traceYdAdjYd - mHd2*traceYeAdjYe)*Sqr(QHd) + traceconjTYeTpTYe*Sqr(Ql) -
      MassU*traceconjTYeTpYe*Sqr(Ql) + traceme2YeAdjYe*Sqr(Ql) +
      traceml2AdjYeYe*Sqr(Ql) + mHd2*traceYeAdjYe*Sqr(Ql) + 3*traceconjTYdTpTYd
      *Sqr(Qq) - 3*MassU*traceconjTYdTpYd*Sqr(Qq) + 3*tracemd2YdAdjYd*Sqr(Qq) +
      3*tracemq2AdjYdYd*Sqr(Qq) + 3*mHd2*traceYdAdjYd*Sqr(Qq))))));
   const double beta_mHd2_2 = Re(0.2*twoLoop*(480*Power(gp,4)*Power(QHd,4
      )*AbsSqr(MassU) + 165*Power(g2,4)*AbsSqr(MassWB) - 30*traceconjTYuTpTYu*
      AbsSqr(Lambdax) - 10*traceconjTYvTpTYv*AbsSqr(Lambdax) - 30*
      tracemq2AdjYuYu*AbsSqr(Lambdax) - 30*tracemu2YuAdjYu*AbsSqr(Lambdax) - 30
      *mHd2*traceYuAdjYu*AbsSqr(Lambdax) - 60*mHu2*traceYuAdjYu*AbsSqr(Lambdax)
      - 30*ms2*traceYuAdjYu*AbsSqr(Lambdax) - 10*mHd2*traceYvAdjYv*AbsSqr(
      Lambdax) - 20*mHu2*traceYvAdjYv*AbsSqr(Lambdax) - 10*ms2*traceYvAdjYv*
      AbsSqr(Lambdax) - 10*traceYvAdjYvconjml2*AbsSqr(Lambdax) - 10*
      traceYvconjmvR2AdjYv*AbsSqr(Lambdax) - 30*traceYuAdjYu*AbsSqr(TLambdax) -
      10*traceYvAdjYv*AbsSqr(TLambdax) - 120*AbsSqr(Lambdax)*AbsSqr(TLambdax)
      - 30*traceAdjYuTYu*Conj(TLambdax)*Lambdax - 10*traceAdjYvTYv*Conj(
      TLambdax)*Lambdax + 18*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 9*MassB*Conj(
      MassWB)*Sqr(g1)*Sqr(g2) - 160*(traceAdjYdTYd - 2*MassG*traceYdAdjYd)*Conj
      (MassG)*Sqr(g3) - 72*Qd*QHd*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 72*Qe*QHd*
      AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 24*QHd*QHu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)
      + 72*QHd*Ql*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 72*QHd*Qq*AbsSqr(MassU)*Sqr(
      g1)*Sqr(gp) + 144*QHd*Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 72*QHd*Qv*AbsSqr
      (MassU)*Sqr(g1)*Sqr(gp) - 36*MassB*Qd*QHd*Conj(MassU)*Sqr(g1)*Sqr(gp) -
      36*MassB*Qe*QHd*Conj(MassU)*Sqr(g1)*Sqr(gp) - 12*MassB*QHd*QHu*Conj(MassU
      )*Sqr(g1)*Sqr(gp) + 36*MassB*QHd*Ql*Conj(MassU)*Sqr(g1)*Sqr(gp) - 36*
      MassB*QHd*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp) + 72*MassB*QHd*Qu*Conj(MassU)*
      Sqr(g1)*Sqr(gp) - 36*MassB*QHd*Qv*Conj(MassU)*Sqr(g1)*Sqr(gp) - 4*Conj(
      MassB)*Sqr(g1)*(-traceAdjYdTYd + 3*traceAdjYeTYe + 2*MassB*traceYdAdjYd -
      6*MassB*traceYeAdjYe + 9*MassU*QHd*Qv*Sqr(gp)) + 120*traceYdAdjYd*AbsSqr
      (MassU)*Sqr(gp)*Sqr(Qd) - 60*traceAdjYdTYd*Conj(MassU)*Sqr(gp)*Sqr(Qd) +
      40*traceYeAdjYe*AbsSqr(MassU)*Sqr(gp)*Sqr(Qe) - 20*traceAdjYeTYe*Conj(
      MassU)*Sqr(gp)*Sqr(Qe) - 120*traceYdAdjYd*AbsSqr(MassU)*Sqr(gp)*Sqr(QHd)
      - 40*traceYeAdjYe*AbsSqr(MassU)*Sqr(gp)*Sqr(QHd) - 20*mHd2*AbsSqr(Lambdax
      )*Sqr(gp)*Sqr(QHd) - 20*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) - 20*ms2*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) - 40*AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(
      gp)*Sqr(QHd) - 20*AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHd) + 60*traceAdjYdTYd*
      Conj(MassU)*Sqr(gp)*Sqr(QHd) + 20*traceAdjYeTYe*Conj(MassU)*Sqr(gp)*Sqr(
      QHd) + 20*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(QHd) + 48*AbsSqr(MassU
      )*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 24*MassB*Conj(MassU)*Sqr(g1)*Sqr(gp)*Sqr(QHd
      ) + 120*AbsSqr(MassU)*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 120*AbsSqr(MassWB)*Sqr(
      g2)*Sqr(gp)*Sqr(QHd) + 60*MassWB*Conj(MassU)*Sqr(g2)*Sqr(gp)*Sqr(QHd) +
      60*MassU*Conj(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(QHd) + 1080*Power(gp,4)*AbsSqr(
      MassU)*Sqr(Qd)*Sqr(QHd) + 360*Power(gp,4)*AbsSqr(MassU)*Sqr(Qe)*Sqr(QHd)
      + 20*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 20*mHu2*AbsSqr(Lambdax)*Sqr(
      gp)*Sqr(QHu) + 20*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 40*AbsSqr(MassU)
      *AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) + 20*AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHu)
      - 20*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(QHu) + 240*Power(gp,4)*
      AbsSqr(MassU)*Sqr(QHd)*Sqr(QHu) + 40*traceYeAdjYe*AbsSqr(MassU)*Sqr(gp)*
      Sqr(Ql) - 20*traceAdjYeTYe*Conj(MassU)*Sqr(gp)*Sqr(Ql) + 720*Power(gp,4)*
      AbsSqr(MassU)*Sqr(QHd)*Sqr(Ql) + 120*traceYdAdjYd*AbsSqr(MassU)*Sqr(gp)*
      Sqr(Qq) - 60*traceAdjYdTYd*Conj(MassU)*Sqr(gp)*Sqr(Qq) + 2160*Power(gp,4)
      *AbsSqr(MassU)*Sqr(QHd)*Sqr(Qq) + 20*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs)
      + 20*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*ms2*AbsSqr(Lambdax)*Sqr(
      gp)*Sqr(Qs) + 40*AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*
      AbsSqr(TLambdax)*Sqr(gp)*Sqr(Qs) - 20*MassU*Conj(TLambdax)*Lambdax*Sqr(gp
      )*Sqr(Qs) + 120*Power(gp,4)*AbsSqr(MassU)*Sqr(QHd)*Sqr(Qs) + 1080*Power(
      gp,4)*AbsSqr(MassU)*Sqr(QHd)*Sqr(Qu) + 360*Power(gp,4)*AbsSqr(MassU)*Sqr(
      QHd)*Sqr(Qv) - 60*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 60*mHu2*Sqr(Conj
      (Lambdax))*Sqr(Lambdax) - 60*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 30*
      traceconjTYuTpYu*Conj(Lambdax)*TLambdax - 10*traceconjTYvTpYv*Conj(
      Lambdax)*TLambdax + 20*Conj(MassU)*Conj(Lambdax)*Sqr(gp)*Sqr(QHd)*
      TLambdax - 20*Conj(MassU)*Conj(Lambdax)*Sqr(gp)*Sqr(QHu)*TLambdax - 20*
      Conj(MassU)*Conj(Lambdax)*Sqr(gp)*Sqr(Qs)*TLambdax));

   beta_mHd2 = beta_mHd2_1 + beta_mHd2_2;


   return beta_mHd2;
}

/**
 * Calculates the three-loop beta function of mHd2.
 *
 * @return three-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

} // namespace flexiblesusy
