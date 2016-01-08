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

// File generated at Fri 8 Jan 2016 15:12:26

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
 * Calculates the one-loop beta function of mHu2.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_one_loop(const Soft_traces& soft_traces) const
{
   const auto QHu = INPUT(QHu);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHu2;

   beta_mHu2 = Re(oneOver16PiSqr*(0.7745966692414834*g1*Tr11 + 2*gp*QHu*
      Tr14 + 6*traceconjTYuTpTYu + 2*traceconjTYvTpTYv + 6*tracemq2AdjYuYu + 6*
      tracemu2YuAdjYu + 6*mHu2*traceYuAdjYu + 2*mHu2*traceYvAdjYv + 2*
      traceYvAdjYvconjml2 + 2*traceYvconjmvR2AdjYv + 2*mHd2*AbsSqr(Lambdax) + 2
      *mHu2*AbsSqr(Lambdax) + 2*ms2*AbsSqr(Lambdax) + 2*AbsSqr(TLambdax) - 1.2*
      AbsSqr(MassB)*Sqr(g1) - 6*AbsSqr(MassWB)*Sqr(g2) - 8*AbsSqr(MassU)*Sqr(gp
      )*Sqr(QHu)));


   return beta_mHu2;
}

/**
 * Calculates the two-loop beta function of mHu2.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_two_loop(const Soft_traces& soft_traces) const
{
   const auto QHu = INPUT(QHu);
   const auto QHd = INPUT(QHd);
   const auto Qs = INPUT(Qs);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Ql = INPUT(Ql);
   const auto Qv = INPUT(Qv);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
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
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeconjTYvTpTYvAdjYe =
      TRACE_STRUCT.traceYeconjTYvTpTYvAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu =
      TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu =
      TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTYvAdjTYv =
      TRACE_STRUCT.traceYvAdjYvTYvAdjTYv;
   const double traceYvAdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYe;
   const double traceYvAdjYvTpTYeconjTYe =
      TRACE_STRUCT.traceYvAdjYvTpTYeconjTYe;
   const double traceYvAdjTYvTYvAdjYv =
      TRACE_STRUCT.traceYvAdjTYvTYvAdjYv;
   const double traceAdjYeTYeconjTYvTpYv =
      TRACE_STRUCT.traceAdjYeTYeconjTYvTpYv;
   const double traceAdjYvTpYeconjTYeTYv =
      TRACE_STRUCT.traceAdjYvTpYeconjTYeTYv;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYvconjml2 =
      TRACE_STRUCT.traceYvAdjYvYvAdjYvconjml2;
   const double traceYvAdjYvYvconjmvR2AdjYv =
      TRACE_STRUCT.traceYvAdjYvYvconjmvR2AdjYv;
   const double traceYvAdjYvconjml2YvAdjYv =
      TRACE_STRUCT.traceYvAdjYvconjml2YvAdjYv;
   const double traceYvAdjYvconjml2TpYeconjYe =
      TRACE_STRUCT.traceYvAdjYvconjml2TpYeconjYe;
   const double traceYvAdjYvTpYeconjme2conjYe =
      TRACE_STRUCT.traceYvAdjYvTpYeconjme2conjYe;
   const double traceYvAdjYvTpYeconjYeconjml2 =
      TRACE_STRUCT.traceYvAdjYvTpYeconjYeconjml2;
   const double traceYvconjmvR2AdjYvYvAdjYv =
      TRACE_STRUCT.traceYvconjmvR2AdjYvYvAdjYv;
   const double traceYvconjmvR2AdjYvTpYeconjYe =
      TRACE_STRUCT.traceYvconjmvR2AdjYvTpYeconjYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHu2;

   const double beta_mHu2_1 = Re(0.04*twoLoop*(3*Conj(MassB)*Sqr(g1)*(261
      *MassB*Sqr(g1) + 5*(3*(2*MassB + MassWB)*Sqr(g2) + 4*QHu*(MassU*(3*Qd + 3
      *Qe - QHd + 2*QHu - 3*Ql + 3*Qq) + 2*MassB*(3*Qd + 3*Qe - QHd + 2*QHu - 3
      *Ql + 3*Qq - 6*Qu))*Sqr(gp))) + 10*(15*Power(g2,4)*Tr22 +
      7.745966692414834*g1*(gp*QHu*(Tr2U114 + Tr2U141) + Tr31) + (3*Tr2U111 + 4
      *traceconjTYuTpTYu - 4*MassB*traceconjTYuTpYu + 6*traceconjTYvTpTYv - 6*
      MassB*traceconjTYvTpYv + 4*tracemq2AdjYuYu + 4*tracemu2YuAdjYu + 4*mHu2*
      traceYuAdjYu + 6*mHu2*traceYvAdjYv + 6*traceYvAdjYvconjml2 + 6*
      traceYvconjmvR2AdjYv)*Sqr(g1) - 5*(-4*gp*QHu*Tr34 +
      traceAdjYeTYeconjTYvTpYv + traceAdjYvTpYeconjTYeTYv + 3*
      tracemd2YdAdjYuYuAdjYd + 3*tracemq2AdjYdYdAdjYuYu + 3*
      tracemq2AdjYuYuAdjYdYd + 18*tracemq2AdjYuYuAdjYuYu + 3*
      tracemu2YuAdjYdYdAdjYu + 18*tracemu2YuAdjYuYuAdjYu + 3*
      traceYdAdjTYuTYuAdjYd + 3*traceYdAdjYuTYuAdjTYd + 3*mHd2*
      traceYdAdjYuYuAdjYd + 3*mHu2*traceYdAdjYuYuAdjYd +
      traceYeconjTYvTpTYvAdjYe + 3*traceYuAdjTYdTYdAdjYu + 18*
      traceYuAdjTYuTYuAdjYu + 3*traceYuAdjYdTYdAdjTYu + 18*
      traceYuAdjYuTYuAdjTYu + 18*mHu2*traceYuAdjYuYuAdjYu + 6*
      traceYvAdjTYvTYvAdjYv + traceYvAdjYvconjml2TpYeconjYe + 3*
      traceYvAdjYvconjml2YvAdjYv + traceYvAdjYvTpTYeconjTYe +
      traceYvAdjYvTpYeconjme2conjYe + mHd2*traceYvAdjYvTpYeconjYe + mHu2*
      traceYvAdjYvTpYeconjYe + traceYvAdjYvTpYeconjYeconjml2 + 6*
      traceYvAdjYvTYvAdjTYv + 6*mHu2*traceYvAdjYvYvAdjYv + 3*
      traceYvAdjYvYvAdjYvconjml2 + 3*traceYvAdjYvYvconjmvR2AdjYv +
      traceYvconjmvR2AdjYvTpYeconjYe + 3*traceYvconjmvR2AdjYvYvAdjYv - 16*
      traceconjTYuTpTYu*Sqr(g3) + 16*MassG*traceconjTYuTpYu*Sqr(g3) - 16*
      tracemq2AdjYuYu*Sqr(g3) - 16*tracemu2YuAdjYu*Sqr(g3) - 16*mHu2*
      traceYuAdjYu*Sqr(g3) - 2*Sqr(gp)*((2*Tr2U144 - 3*traceconjTYuTpTYu + 3*
      MassU*traceconjTYuTpYu - traceconjTYvTpTYv + MassU*traceconjTYvTpYv - 3*
      tracemq2AdjYuYu - 3*tracemu2YuAdjYu - 3*mHu2*traceYuAdjYu - mHu2*
      traceYvAdjYv - traceYvAdjYvconjml2 - traceYvconjmvR2AdjYv)*Sqr(QHu) +
      traceconjTYvTpTYv*Sqr(Ql) - MassU*traceconjTYvTpYv*Sqr(Ql) + mHu2*
      traceYvAdjYv*Sqr(Ql) + traceYvAdjYvconjml2*Sqr(Ql) + traceYvconjmvR2AdjYv
      *Sqr(Ql) + 3*(traceconjTYuTpTYu - MassU*traceconjTYuTpYu +
      tracemq2AdjYuYu + tracemu2YuAdjYu + mHu2*traceYuAdjYu)*Sqr(Qq) + 3*
      traceconjTYuTpTYu*Sqr(Qu) - 3*MassU*traceconjTYuTpYu*Sqr(Qu) + 3*
      tracemq2AdjYuYu*Sqr(Qu) + 3*tracemu2YuAdjYu*Sqr(Qu) + 3*mHu2*traceYuAdjYu
      *Sqr(Qu) + traceconjTYvTpTYv*Sqr(Qv) - MassU*traceconjTYvTpYv*Sqr(Qv) +
      mHu2*traceYvAdjYv*Sqr(Qv) + traceYvAdjYvconjml2*Sqr(Qv) +
      traceYvconjmvR2AdjYv*Sqr(Qv))))));
   const double beta_mHu2_2 = Re(0.2*twoLoop*(480*Power(gp,4)*Power(QHu,4
      )*AbsSqr(MassU) + 165*Power(g2,4)*AbsSqr(MassWB) - 30*traceconjTYdTpTYd*
      AbsSqr(Lambdax) - 10*traceconjTYeTpTYe*AbsSqr(Lambdax) - 30*
      tracemd2YdAdjYd*AbsSqr(Lambdax) - 10*traceme2YeAdjYe*AbsSqr(Lambdax) - 10
      *traceml2AdjYeYe*AbsSqr(Lambdax) - 30*tracemq2AdjYdYd*AbsSqr(Lambdax) -
      60*mHd2*traceYdAdjYd*AbsSqr(Lambdax) - 30*mHu2*traceYdAdjYd*AbsSqr(
      Lambdax) - 30*ms2*traceYdAdjYd*AbsSqr(Lambdax) - 20*mHd2*traceYeAdjYe*
      AbsSqr(Lambdax) - 10*mHu2*traceYeAdjYe*AbsSqr(Lambdax) - 10*ms2*
      traceYeAdjYe*AbsSqr(Lambdax) - 30*traceYdAdjYd*AbsSqr(TLambdax) - 10*
      traceYeAdjYe*AbsSqr(TLambdax) - 120*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 30
      *traceAdjYdTYd*Conj(TLambdax)*Lambdax - 10*traceAdjYeTYe*Conj(TLambdax)*
      Lambdax + 18*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 9*MassB*Conj(MassWB)*Sqr(g1
      )*Sqr(g2) - 160*(traceAdjYuTYu - 2*MassG*traceYuAdjYu)*Conj(MassG)*Sqr(g3
      ) + 72*Qd*QHu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) + 72*Qe*QHu*AbsSqr(MassU)*Sqr
      (g1)*Sqr(gp) - 24*QHd*QHu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 72*QHu*Ql*
      AbsSqr(MassU)*Sqr(g1)*Sqr(gp) + 72*QHu*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) -
      144*QHu*Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) + 72*QHu*Qv*AbsSqr(MassU)*Sqr(
      g1)*Sqr(gp) + 36*MassB*Qd*QHu*Conj(MassU)*Sqr(g1)*Sqr(gp) + 36*MassB*Qe*
      QHu*Conj(MassU)*Sqr(g1)*Sqr(gp) - 12*MassB*QHd*QHu*Conj(MassU)*Sqr(g1)*
      Sqr(gp) - 36*MassB*QHu*Ql*Conj(MassU)*Sqr(g1)*Sqr(gp) + 36*MassB*QHu*Qq*
      Conj(MassU)*Sqr(g1)*Sqr(gp) - 72*MassB*QHu*Qu*Conj(MassU)*Sqr(g1)*Sqr(gp)
      + 36*MassB*QHu*Qv*Conj(MassU)*Sqr(g1)*Sqr(gp) - 4*Conj(MassB)*Sqr(g1)*(2
      *traceAdjYuTYu + 3*traceAdjYvTYv - 4*MassB*traceYuAdjYu - 6*MassB*
      traceYvAdjYv + 9*QHu*(2*MassU*Qu - 2*MassB*Qv - MassU*Qv)*Sqr(gp)) + 20*
      mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 20*mHu2*AbsSqr(Lambdax)*Sqr(gp)*
      Sqr(QHd) + 20*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 40*AbsSqr(MassU)*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 20*AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHd) -
      20*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(QHd) - 120*traceYuAdjYu*
      AbsSqr(MassU)*Sqr(gp)*Sqr(QHu) - 40*traceYvAdjYv*AbsSqr(MassU)*Sqr(gp)*
      Sqr(QHu) - 20*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) - 20*mHu2*AbsSqr(
      Lambdax)*Sqr(gp)*Sqr(QHu) - 20*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) - 40*
      AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) - 20*AbsSqr(TLambdax)*Sqr(
      gp)*Sqr(QHu) + 60*traceAdjYuTYu*Conj(MassU)*Sqr(gp)*Sqr(QHu) + 20*
      traceAdjYvTYv*Conj(MassU)*Sqr(gp)*Sqr(QHu) + 20*MassU*Conj(TLambdax)*
      Lambdax*Sqr(gp)*Sqr(QHu) + 48*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 24
      *MassB*Conj(MassU)*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 120*AbsSqr(MassU)*Sqr(g2)*
      Sqr(gp)*Sqr(QHu) + 120*AbsSqr(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 60*
      MassWB*Conj(MassU)*Sqr(g2)*Sqr(gp)*Sqr(QHu) + 60*MassU*Conj(MassWB)*Sqr(
      g2)*Sqr(gp)*Sqr(QHu) + 1080*Power(gp,4)*AbsSqr(MassU)*Sqr(Qd)*Sqr(QHu) +
      360*Power(gp,4)*AbsSqr(MassU)*Sqr(Qe)*Sqr(QHu) + 240*Power(gp,4)*AbsSqr(
      MassU)*Sqr(QHd)*Sqr(QHu) + 40*traceYvAdjYv*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql)
      - 20*traceAdjYvTYv*Conj(MassU)*Sqr(gp)*Sqr(Ql) + 720*Power(gp,4)*AbsSqr(
      MassU)*Sqr(QHu)*Sqr(Ql) + 120*traceYuAdjYu*AbsSqr(MassU)*Sqr(gp)*Sqr(Qq)
      - 60*traceAdjYuTYu*Conj(MassU)*Sqr(gp)*Sqr(Qq) + 2160*Power(gp,4)*AbsSqr(
      MassU)*Sqr(QHu)*Sqr(Qq) + 20*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*
      mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr
      (Qs) + 40*AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*AbsSqr(
      TLambdax)*Sqr(gp)*Sqr(Qs) - 20*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(
      Qs) + 120*Power(gp,4)*AbsSqr(MassU)*Sqr(QHu)*Sqr(Qs) + 120*traceYuAdjYu*
      AbsSqr(MassU)*Sqr(gp)*Sqr(Qu) - 60*traceAdjYuTYu*Conj(MassU)*Sqr(gp)*Sqr(
      Qu) + 1080*Power(gp,4)*AbsSqr(MassU)*Sqr(QHu)*Sqr(Qu) + 40*traceYvAdjYv*
      AbsSqr(MassU)*Sqr(gp)*Sqr(Qv) - 20*traceAdjYvTYv*Conj(MassU)*Sqr(gp)*Sqr(
      Qv) + 360*Power(gp,4)*AbsSqr(MassU)*Sqr(QHu)*Sqr(Qv) - 60*mHd2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) - 60*mHu2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 60*ms2
      *Sqr(Conj(Lambdax))*Sqr(Lambdax) - 30*traceconjTYdTpYd*Conj(Lambdax)*
      TLambdax - 10*traceconjTYeTpYe*Conj(Lambdax)*TLambdax - 20*Conj(MassU)*
      Conj(Lambdax)*Sqr(gp)*Sqr(QHd)*TLambdax + 20*Conj(MassU)*Conj(Lambdax)*
      Sqr(gp)*Sqr(QHu)*TLambdax - 20*Conj(MassU)*Conj(Lambdax)*Sqr(gp)*Sqr(Qs)*
      TLambdax));

   beta_mHu2 = beta_mHu2_1 + beta_mHu2_2;


   return beta_mHu2;
}

/**
 * Calculates the three-loop beta function of mHu2.
 *
 * @return three-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

} // namespace flexiblesusy
