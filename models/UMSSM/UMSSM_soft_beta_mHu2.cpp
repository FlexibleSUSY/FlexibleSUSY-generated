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

// File generated at Fri 10 Apr 2020 20:18:29

#include "UMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mHu2.
 *
 * @return 1-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_1_loop(const Soft_traces& soft_traces) const
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

   beta_mHu2 = Re(0.2*oneOver16PiSqr*(3.872983346207417*g1*Tr11 + 10*gp*QHu*
      Tr14 + 30*traceconjTYuTpTYu + 10*traceconjTYvTpTYv + 30*tracemq2AdjYuYu +
      30*tracemu2YuAdjYu + 30*mHu2*traceYuAdjYu + 10*mHu2*traceYvAdjYv + 10*
      traceYvAdjYvconjml2 + 10*traceYvconjmvR2AdjYv + 10*mHd2*AbsSqr(Lambdax) +
      10*mHu2*AbsSqr(Lambdax) + 10*ms2*AbsSqr(Lambdax) + 10*AbsSqr(TLambdax) -
      6*AbsSqr(MassB)*Sqr(g1) - 30*AbsSqr(MassWB)*Sqr(g2) - 40*AbsSqr(MassU)*
      Sqr(gp)*Sqr(QHu)));


   return beta_mHu2;
}

/**
 * Calculates the 2-loop beta function of mHu2.
 *
 * @return 2-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_2_loop(const Soft_traces& soft_traces) const
{
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Qs = INPUT(Qs);
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
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
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd = TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd = TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeconjTYvTpTYvAdjYe = TRACE_STRUCT.
      traceYeconjTYvTpTYvAdjYe;
   const double traceYuAdjYdTYdAdjTYu = TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu = TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu = TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu = TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double traceYvAdjYvYvAdjYv = TRACE_STRUCT.traceYvAdjYvYvAdjYv;
   const double traceYvAdjYvTYvAdjTYv = TRACE_STRUCT.traceYvAdjYvTYvAdjTYv;
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;
   const double traceYvAdjYvTpTYeconjTYe = TRACE_STRUCT.
      traceYvAdjYvTpTYeconjTYe;
   const double traceYvAdjTYvTYvAdjYv = TRACE_STRUCT.traceYvAdjTYvTYvAdjYv;
   const double traceAdjYeTYeconjTYvTpYv = TRACE_STRUCT.
      traceAdjYeTYeconjTYvTpYv;
   const double traceAdjYvTpYeconjTYeTYv = TRACE_STRUCT.
      traceAdjYvTpYeconjTYeTYv;
   const double tracemd2YdAdjYuYuAdjYd = TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double tracemq2AdjYdYdAdjYuYu = TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd = TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu = TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu = TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu = TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double traceYvAdjYvYvAdjYvconjml2 = TRACE_STRUCT.
      traceYvAdjYvYvAdjYvconjml2;
   const double traceYvAdjYvYvconjmvR2AdjYv = TRACE_STRUCT.
      traceYvAdjYvYvconjmvR2AdjYv;
   const double traceYvAdjYvconjml2YvAdjYv = TRACE_STRUCT.
      traceYvAdjYvconjml2YvAdjYv;
   const double traceYvAdjYvconjml2TpYeconjYe = TRACE_STRUCT.
      traceYvAdjYvconjml2TpYeconjYe;
   const double traceYvAdjYvTpYeconjme2conjYe = TRACE_STRUCT.
      traceYvAdjYvTpYeconjme2conjYe;
   const double traceYvAdjYvTpYeconjYeconjml2 = TRACE_STRUCT.
      traceYvAdjYvTpYeconjYeconjml2;
   const double traceYvconjmvR2AdjYvYvAdjYv = TRACE_STRUCT.
      traceYvconjmvR2AdjYvYvAdjYv;
   const double traceYvconjmvR2AdjYvTpYeconjYe = TRACE_STRUCT.
      traceYvconjmvR2AdjYvTpYeconjYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHu2;

   const double beta_mHu2_1 = Re(0.04*twoLoop*(77.45966692414834*g1*gp*QHu*
      Tr2U114 + 77.45966692414834*g1*gp*QHu*Tr2U141 + 77.45966692414834*g1*Tr31
       + 200*gp*QHu*Tr34 - 50*traceAdjYeTYeconjTYvTpYv - 50*
      traceAdjYvTpYeconjTYeTYv - 150*tracemd2YdAdjYuYuAdjYd - 150*
      tracemq2AdjYdYdAdjYuYu - 150*tracemq2AdjYuYuAdjYdYd - 900*
      tracemq2AdjYuYuAdjYuYu - 150*tracemu2YuAdjYdYdAdjYu - 900*
      tracemu2YuAdjYuYuAdjYu - 150*traceYdAdjTYuTYuAdjYd - 150*
      traceYdAdjYuTYuAdjTYd - 150*mHd2*traceYdAdjYuYuAdjYd - 150*mHu2*
      traceYdAdjYuYuAdjYd - 50*traceYeconjTYvTpTYvAdjYe - 150*
      traceYuAdjTYdTYdAdjYu - 900*traceYuAdjTYuTYuAdjYu - 150*
      traceYuAdjYdTYdAdjTYu - 900*traceYuAdjYuTYuAdjTYu - 900*mHu2*
      traceYuAdjYuYuAdjYu - 300*traceYvAdjTYvTYvAdjYv - 50*
      traceYvAdjYvconjml2TpYeconjYe - 150*traceYvAdjYvconjml2YvAdjYv - 50*
      traceYvAdjYvTpTYeconjTYe - 50*traceYvAdjYvTpYeconjme2conjYe - 50*mHd2*
      traceYvAdjYvTpYeconjYe - 50*mHu2*traceYvAdjYvTpYeconjYe - 50*
      traceYvAdjYvTpYeconjYeconjml2 - 300*traceYvAdjYvTYvAdjTYv - 300*mHu2*
      traceYvAdjYvYvAdjYv - 150*traceYvAdjYvYvAdjYvconjml2 - 150*
      traceYvAdjYvYvconjmvR2AdjYv - 50*traceYvconjmvR2AdjYvTpYeconjYe - 150*
      traceYvconjmvR2AdjYvYvAdjYv + 621*AbsSqr(MassB)*Quad(g1) + 150*Tr22*Quad(
      g2) + 30*Tr2U111*Sqr(g1) + 40*traceconjTYuTpTYu*Sqr(g1) - 40*MassB*
      traceconjTYuTpYu*Sqr(g1) + 40*tracemq2AdjYuYu*Sqr(g1) + 40*
      tracemu2YuAdjYu*Sqr(g1) + 40*mHu2*traceYuAdjYu*Sqr(g1) + 80*traceYuAdjYu*
      AbsSqr(MassB)*Sqr(g1) - 40*traceAdjYuTYu*Conj(MassB)*Sqr(g1) + 90*AbsSqr(
      MassB)*Sqr(g1)*Sqr(g2) + 45*MassWB*Conj(MassB)*Sqr(g1)*Sqr(g2) + 800*
      traceconjTYuTpTYu*Sqr(g3) - 800*MassG*traceconjTYuTpYu*Sqr(g3) + 800*
      tracemq2AdjYuYu*Sqr(g3) + 800*tracemu2YuAdjYu*Sqr(g3) + 800*mHu2*
      traceYuAdjYu*Sqr(g3) + 1600*traceYuAdjYu*AbsSqr(MassG)*Sqr(g3) - 800*
      traceAdjYuTYu*Conj(MassG)*Sqr(g3) + 360*Qd*QHu*AbsSqr(MassB)*Sqr(g1)*Sqr(
      gp) + 360*Qe*QHu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) - 120*QHd*QHu*AbsSqr(MassB
      )*Sqr(g1)*Sqr(gp) - 360*QHu*Ql*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 360*QHu*Qq
      *AbsSqr(MassB)*Sqr(g1)*Sqr(gp) - 720*QHu*Qu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)
      + 180*MassU*Qd*QHu*Conj(MassB)*Sqr(g1)*Sqr(gp) + 180*MassU*Qe*QHu*Conj(
      MassB)*Sqr(g1)*Sqr(gp) - 60*MassU*QHd*QHu*Conj(MassB)*Sqr(g1)*Sqr(gp) -
      180*MassU*QHu*Ql*Conj(MassB)*Sqr(g1)*Sqr(gp) + 180*MassU*QHu*Qq*Conj(
      MassB)*Sqr(g1)*Sqr(gp) - 360*MassU*QHu*Qu*Conj(MassB)*Sqr(g1)*Sqr(gp) +
      200*Tr2U144*Sqr(gp)*Sqr(QHu) - 300*traceconjTYuTpTYu*Sqr(gp)*Sqr(QHu) +
      300*MassU*traceconjTYuTpYu*Sqr(gp)*Sqr(QHu) - 100*traceconjTYvTpTYv*Sqr(
      gp)*Sqr(QHu) + 100*MassU*traceconjTYvTpYv*Sqr(gp)*Sqr(QHu) - 300*
      tracemq2AdjYuYu*Sqr(gp)*Sqr(QHu) - 300*tracemu2YuAdjYu*Sqr(gp)*Sqr(QHu) -
      300*mHu2*traceYuAdjYu*Sqr(gp)*Sqr(QHu) - 100*mHu2*traceYvAdjYv*Sqr(gp)*
      Sqr(QHu) - 100*traceYvAdjYvconjml2*Sqr(gp)*Sqr(QHu) - 100*
      traceYvconjmvR2AdjYv*Sqr(gp)*Sqr(QHu) + 240*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)
      *Sqr(QHu) + 120*MassU*Conj(MassB)*Sqr(g1)*Sqr(gp)*Sqr(QHu) + 100*
      traceconjTYvTpTYv*Sqr(gp)*Sqr(Ql) - 100*MassU*traceconjTYvTpYv*Sqr(gp)*
      Sqr(Ql) + 100*mHu2*traceYvAdjYv*Sqr(gp)*Sqr(Ql) + 100*traceYvAdjYvconjml2
      *Sqr(gp)*Sqr(Ql) + 100*traceYvconjmvR2AdjYv*Sqr(gp)*Sqr(Ql) + 300*
      traceconjTYuTpTYu*Sqr(gp)*Sqr(Qq) - 300*MassU*traceconjTYuTpYu*Sqr(gp)*
      Sqr(Qq) + 300*tracemq2AdjYuYu*Sqr(gp)*Sqr(Qq) + 300*tracemu2YuAdjYu*Sqr(
      gp)*Sqr(Qq) + 300*mHu2*traceYuAdjYu*Sqr(gp)*Sqr(Qq) + 300*
      traceconjTYuTpTYu*Sqr(gp)*Sqr(Qu) - 300*MassU*traceconjTYuTpYu*Sqr(gp)*
      Sqr(Qu) + 300*tracemq2AdjYuYu*Sqr(gp)*Sqr(Qu) + 300*tracemu2YuAdjYu*Sqr(
      gp)*Sqr(Qu) + 300*mHu2*traceYuAdjYu*Sqr(gp)*Sqr(Qu) + 100*
      traceconjTYvTpTYv*Sqr(gp)*Sqr(Qv) - 100*MassU*traceconjTYvTpYv*Sqr(gp)*
      Sqr(Qv) + 100*mHu2*traceYvAdjYv*Sqr(gp)*Sqr(Qv) + 100*traceYvAdjYvconjml2
      *Sqr(gp)*Sqr(Qv) + 100*traceYvconjmvR2AdjYv*Sqr(gp)*Sqr(Qv)));
   const double beta_mHu2_2 = Re(0.2*twoLoop*(-30*traceconjTYdTpTYd*AbsSqr(
      Lambdax) - 10*traceconjTYeTpTYe*AbsSqr(Lambdax) - 30*tracemd2YdAdjYd*
      AbsSqr(Lambdax) - 10*traceme2YeAdjYe*AbsSqr(Lambdax) - 10*traceml2AdjYeYe
      *AbsSqr(Lambdax) - 30*tracemq2AdjYdYd*AbsSqr(Lambdax) - 60*mHd2*
      traceYdAdjYd*AbsSqr(Lambdax) - 30*mHu2*traceYdAdjYd*AbsSqr(Lambdax) - 30*
      ms2*traceYdAdjYd*AbsSqr(Lambdax) - 20*mHd2*traceYeAdjYe*AbsSqr(Lambdax) -
      10*mHu2*traceYeAdjYe*AbsSqr(Lambdax) - 10*ms2*traceYeAdjYe*AbsSqr(Lambdax
      ) - 30*traceYdAdjYd*AbsSqr(TLambdax) - 10*traceYeAdjYe*AbsSqr(TLambdax) -
      120*AbsSqr(Lambdax)*AbsSqr(TLambdax) - 30*traceAdjYdTYd*Conj(TLambdax)*
      Lambdax - 10*traceAdjYeTYe*Conj(TLambdax)*Lambdax + 165*AbsSqr(MassWB)*
      Quad(g2) + 480*AbsSqr(MassU)*Quad(gp)*Quad(QHu) + 18*AbsSqr(MassWB)*Sqr(
      g1)*Sqr(g2) + 9*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 72*Qd*QHu*AbsSqr(
      MassU)*Sqr(g1)*Sqr(gp) + 72*Qe*QHu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 24*QHd
      *QHu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 72*QHu*Ql*AbsSqr(MassU)*Sqr(g1)*Sqr(
      gp) + 72*QHu*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 144*QHu*Qu*AbsSqr(MassU)*
      Sqr(g1)*Sqr(gp) + 36*MassB*Qd*QHu*Conj(MassU)*Sqr(g1)*Sqr(gp) + 36*MassB*
      Qe*QHu*Conj(MassU)*Sqr(g1)*Sqr(gp) - 12*MassB*QHd*QHu*Conj(MassU)*Sqr(g1)
      *Sqr(gp) - 36*MassB*QHu*Ql*Conj(MassU)*Sqr(g1)*Sqr(gp) + 36*MassB*QHu*Qq*
      Conj(MassU)*Sqr(g1)*Sqr(gp) - 72*MassB*QHu*Qu*Conj(MassU)*Sqr(g1)*Sqr(gp)
      + 20*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 20*mHu2*AbsSqr(Lambdax)*Sqr(
      gp)*Sqr(QHd) + 20*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 40*AbsSqr(MassU)
      *AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 20*AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHd)
      - 20*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(QHd) - 120*traceYuAdjYu*
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
      g2)*Sqr(gp)*Sqr(QHu) + 1080*AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr(QHu) + 360
      *AbsSqr(MassU)*Quad(gp)*Sqr(Qe)*Sqr(QHu) + 240*AbsSqr(MassU)*Quad(gp)*Sqr
      (QHd)*Sqr(QHu) + 40*traceYvAdjYv*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql) - 20*
      traceAdjYvTYv*Conj(MassU)*Sqr(gp)*Sqr(Ql) + 720*AbsSqr(MassU)*Quad(gp)*
      Sqr(QHu)*Sqr(Ql) + 120*traceYuAdjYu*AbsSqr(MassU)*Sqr(gp)*Sqr(Qq) - 60*
      traceAdjYuTYu*Conj(MassU)*Sqr(gp)*Sqr(Qq) + 2160*AbsSqr(MassU)*Quad(gp)*
      Sqr(QHu)*Sqr(Qq) + 20*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*mHu2*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs)
      + 40*AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) + 20*AbsSqr(TLambdax)*
      Sqr(gp)*Sqr(Qs) - 20*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(Qs) + 120*
      AbsSqr(MassU)*Quad(gp)*Sqr(QHu)*Sqr(Qs) + 120*traceYuAdjYu*AbsSqr(MassU)*
      Sqr(gp)*Sqr(Qu) - 60*traceAdjYuTYu*Conj(MassU)*Sqr(gp)*Sqr(Qu) + 1080*
      AbsSqr(MassU)*Quad(gp)*Sqr(QHu)*Sqr(Qu) + 40*traceYvAdjYv*AbsSqr(MassU)*
      Sqr(gp)*Sqr(Qv) - 20*traceAdjYvTYv*Conj(MassU)*Sqr(gp)*Sqr(Qv) + 360*
      AbsSqr(MassU)*Quad(gp)*Sqr(QHu)*Sqr(Qv) - 60*mHd2*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) - 60*mHu2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 60*ms2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) - 30*traceconjTYdTpYd*Conj(Lambdax)*TLambdax - 10*
      traceconjTYeTpYe*Conj(Lambdax)*TLambdax - 20*Conj(MassU)*Conj(Lambdax)*
      Sqr(gp)*Sqr(QHd)*TLambdax + 20*Conj(MassU)*Conj(Lambdax)*Sqr(gp)*Sqr(QHu)
      *TLambdax - 20*Conj(MassU)*Conj(Lambdax)*Sqr(gp)*Sqr(Qs)*TLambdax));

   beta_mHu2 = beta_mHu2_1 + beta_mHu2_2;


   return beta_mHu2;
}

/**
 * Calculates the 3-loop beta function of mHu2.
 *
 * @return 3-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

/**
 * Calculates the 4-loop beta function of mHu2.
 *
 * @return 4-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

/**
 * Calculates the 5-loop beta function of mHu2.
 *
 * @return 5-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHu2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

} // namespace flexiblesusy
