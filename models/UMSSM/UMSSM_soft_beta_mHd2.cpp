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

// File generated at Wed 16 Oct 2019 22:27:15

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
 * Calculates the 1-loop beta function of mHd2.
 *
 * @return 1-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_1_loop(const Soft_traces& soft_traces) const
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

   beta_mHd2 = Re(0.2*oneOver16PiSqr*(-3.872983346207417*g1*Tr11 + 10*gp*QHd*
      Tr14 + 30*traceconjTYdTpTYd + 10*traceconjTYeTpTYe + 30*tracemd2YdAdjYd +
      10*traceme2YeAdjYe + 10*traceml2AdjYeYe + 30*tracemq2AdjYdYd + 30*mHd2*
      traceYdAdjYd + 10*mHd2*traceYeAdjYe + 10*mHd2*AbsSqr(Lambdax) + 10*mHu2*
      AbsSqr(Lambdax) + 10*ms2*AbsSqr(Lambdax) + 10*AbsSqr(TLambdax) - 6*AbsSqr
      (MassB)*Sqr(g1) - 30*AbsSqr(MassWB)*Sqr(g2) - 40*AbsSqr(MassU)*Sqr(gp)*
      Sqr(QHd)));


   return beta_mHd2;
}

/**
 * Calculates the 2-loop beta function of mHd2.
 *
 * @return 2-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_2_loop(const Soft_traces& soft_traces) const
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd = TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd = TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd = TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd = TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe = TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe = TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYeconjTYvTpTYvAdjYe = TRACE_STRUCT.
      traceYeconjTYvTpTYvAdjYe;
   const double traceYuAdjYdTYdAdjTYu = TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu = TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYvAdjYvTpYeconjYe = TRACE_STRUCT.traceYvAdjYvTpYeconjYe;
   const double traceYvAdjYvTpTYeconjTYe = TRACE_STRUCT.
      traceYvAdjYvTpTYeconjTYe;
   const double traceAdjYeTYeconjTYvTpYv = TRACE_STRUCT.
      traceAdjYeTYeconjTYvTpYv;
   const double traceAdjYvTpYeconjTYeTYv = TRACE_STRUCT.
      traceAdjYvTpYeconjTYeTYv;
   const double tracemd2YdAdjYdYdAdjYd = TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd = TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe = TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe = TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd = TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu = TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd = TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu = TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double traceYvAdjYvconjml2TpYeconjYe = TRACE_STRUCT.
      traceYvAdjYvconjml2TpYeconjYe;
   const double traceYvAdjYvTpYeconjme2conjYe = TRACE_STRUCT.
      traceYvAdjYvTpYeconjme2conjYe;
   const double traceYvAdjYvTpYeconjYeconjml2 = TRACE_STRUCT.
      traceYvAdjYvTpYeconjYeconjml2;
   const double traceYvconjmvR2AdjYvTpYeconjYe = TRACE_STRUCT.
      traceYvconjmvR2AdjYvTpYeconjYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHd2;

   const double beta_mHd2_1 = Re(0.04*twoLoop*(-77.45966692414834*g1*gp*QHd*
      Tr2U114 - 77.45966692414834*g1*gp*QHd*Tr2U141 - 77.45966692414834*g1*Tr31
       + 200*gp*QHd*Tr34 - 50*traceAdjYeTYeconjTYvTpYv - 50*
      traceAdjYvTpYeconjTYeTYv - 900*tracemd2YdAdjYdYdAdjYd - 150*
      tracemd2YdAdjYuYuAdjYd - 300*traceme2YeAdjYeYeAdjYe - 300*
      traceml2AdjYeYeAdjYeYe - 900*tracemq2AdjYdYdAdjYdYd - 150*
      tracemq2AdjYdYdAdjYuYu - 150*tracemq2AdjYuYuAdjYdYd - 150*
      tracemu2YuAdjYdYdAdjYu - 900*traceYdAdjTYdTYdAdjYd - 150*
      traceYdAdjTYuTYuAdjYd - 900*traceYdAdjYdTYdAdjTYd - 900*mHd2*
      traceYdAdjYdYdAdjYd - 150*traceYdAdjYuTYuAdjTYd - 150*mHd2*
      traceYdAdjYuYuAdjYd - 150*mHu2*traceYdAdjYuYuAdjYd - 300*
      traceYeAdjTYeTYeAdjYe - 300*traceYeAdjYeTYeAdjTYe - 300*mHd2*
      traceYeAdjYeYeAdjYe - 50*traceYeconjTYvTpTYvAdjYe - 150*
      traceYuAdjTYdTYdAdjYu - 150*traceYuAdjYdTYdAdjTYu - 50*
      traceYvAdjYvconjml2TpYeconjYe - 50*traceYvAdjYvTpTYeconjTYe - 50*
      traceYvAdjYvTpYeconjme2conjYe - 50*mHd2*traceYvAdjYvTpYeconjYe - 50*mHu2*
      traceYvAdjYvTpYeconjYe - 50*traceYvAdjYvTpYeconjYeconjml2 - 50*
      traceYvconjmvR2AdjYvTpYeconjYe + 621*AbsSqr(MassB)*Quad(g1) + 150*Tr22*
      Quad(g2) + 30*Tr2U111*Sqr(g1) - 20*traceconjTYdTpTYd*Sqr(g1) + 20*MassB*
      traceconjTYdTpYd*Sqr(g1) + 60*traceconjTYeTpTYe*Sqr(g1) - 60*MassB*
      traceconjTYeTpYe*Sqr(g1) - 20*tracemd2YdAdjYd*Sqr(g1) + 60*
      traceme2YeAdjYe*Sqr(g1) + 60*traceml2AdjYeYe*Sqr(g1) - 20*tracemq2AdjYdYd
      *Sqr(g1) - 20*mHd2*traceYdAdjYd*Sqr(g1) + 60*mHd2*traceYeAdjYe*Sqr(g1) +
      20*traceAdjYdTYd*Conj(MassB)*Sqr(g1) + 90*AbsSqr(MassB)*Sqr(g1)*Sqr(g2) +
      45*MassWB*Conj(MassB)*Sqr(g1)*Sqr(g2) + 800*traceconjTYdTpTYd*Sqr(g3) -
      800*MassG*traceconjTYdTpYd*Sqr(g3) + 800*tracemd2YdAdjYd*Sqr(g3) + 800*
      tracemq2AdjYdYd*Sqr(g3) + 800*mHd2*traceYdAdjYd*Sqr(g3) - 360*Qd*QHd*
      AbsSqr(MassB)*Sqr(g1)*Sqr(gp) - 360*Qe*QHd*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)
      - 120*QHd*QHu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 360*QHd*Ql*AbsSqr(MassB)*
      Sqr(g1)*Sqr(gp) - 360*QHd*Qq*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 720*QHd*Qu*
      AbsSqr(MassB)*Sqr(g1)*Sqr(gp) - 180*MassU*Qd*QHd*Conj(MassB)*Sqr(g1)*Sqr(
      gp) - 180*MassU*Qe*QHd*Conj(MassB)*Sqr(g1)*Sqr(gp) - 60*MassU*QHd*QHu*
      Conj(MassB)*Sqr(g1)*Sqr(gp) + 180*MassU*QHd*Ql*Conj(MassB)*Sqr(g1)*Sqr(gp
      ) - 180*MassU*QHd*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp) + 360*MassU*QHd*Qu*Conj(
      MassB)*Sqr(g1)*Sqr(gp) + 300*traceconjTYdTpTYd*Sqr(gp)*Sqr(Qd) - 300*
      MassU*traceconjTYdTpYd*Sqr(gp)*Sqr(Qd) + 300*tracemd2YdAdjYd*Sqr(gp)*Sqr(
      Qd) + 300*tracemq2AdjYdYd*Sqr(gp)*Sqr(Qd) + 300*mHd2*traceYdAdjYd*Sqr(gp)
      *Sqr(Qd) + 100*traceconjTYeTpTYe*Sqr(gp)*Sqr(Qe) - 100*MassU*
      traceconjTYeTpYe*Sqr(gp)*Sqr(Qe) + 100*traceme2YeAdjYe*Sqr(gp)*Sqr(Qe) +
      100*traceml2AdjYeYe*Sqr(gp)*Sqr(Qe) + 100*mHd2*traceYeAdjYe*Sqr(gp)*Sqr(
      Qe) + 200*Tr2U144*Sqr(gp)*Sqr(QHd) - 300*traceconjTYdTpTYd*Sqr(gp)*Sqr(
      QHd) + 300*MassU*traceconjTYdTpYd*Sqr(gp)*Sqr(QHd) - 100*
      traceconjTYeTpTYe*Sqr(gp)*Sqr(QHd) + 100*MassU*traceconjTYeTpYe*Sqr(gp)*
      Sqr(QHd) - 300*tracemd2YdAdjYd*Sqr(gp)*Sqr(QHd) - 100*traceme2YeAdjYe*Sqr
      (gp)*Sqr(QHd) - 100*traceml2AdjYeYe*Sqr(gp)*Sqr(QHd) - 300*
      tracemq2AdjYdYd*Sqr(gp)*Sqr(QHd) - 300*mHd2*traceYdAdjYd*Sqr(gp)*Sqr(QHd)
      - 100*mHd2*traceYeAdjYe*Sqr(gp)*Sqr(QHd) + 240*AbsSqr(MassB)*Sqr(g1)*Sqr(
      gp)*Sqr(QHd) + 120*MassU*Conj(MassB)*Sqr(g1)*Sqr(gp)*Sqr(QHd) + 100*
      traceconjTYeTpTYe*Sqr(gp)*Sqr(Ql) - 100*MassU*traceconjTYeTpYe*Sqr(gp)*
      Sqr(Ql) + 100*traceme2YeAdjYe*Sqr(gp)*Sqr(Ql) + 100*traceml2AdjYeYe*Sqr(
      gp)*Sqr(Ql) + 100*mHd2*traceYeAdjYe*Sqr(gp)*Sqr(Ql) + 300*
      traceconjTYdTpTYd*Sqr(gp)*Sqr(Qq) - 300*MassU*traceconjTYdTpYd*Sqr(gp)*
      Sqr(Qq) + 300*tracemd2YdAdjYd*Sqr(gp)*Sqr(Qq) + 300*tracemq2AdjYdYd*Sqr(
      gp)*Sqr(Qq) + 300*mHd2*traceYdAdjYd*Sqr(gp)*Sqr(Qq)));
   const double beta_mHd2_2 = Re(-0.2*twoLoop*(30*traceconjTYuTpTYu*AbsSqr(
      Lambdax) + 10*traceconjTYvTpTYv*AbsSqr(Lambdax) + 30*tracemq2AdjYuYu*
      AbsSqr(Lambdax) + 30*tracemu2YuAdjYu*AbsSqr(Lambdax) + 30*mHd2*
      traceYuAdjYu*AbsSqr(Lambdax) + 60*mHu2*traceYuAdjYu*AbsSqr(Lambdax) + 30*
      ms2*traceYuAdjYu*AbsSqr(Lambdax) + 10*mHd2*traceYvAdjYv*AbsSqr(Lambdax) +
      20*mHu2*traceYvAdjYv*AbsSqr(Lambdax) + 10*ms2*traceYvAdjYv*AbsSqr(Lambdax
      ) + 10*traceYvAdjYvconjml2*AbsSqr(Lambdax) + 10*traceYvconjmvR2AdjYv*
      AbsSqr(Lambdax) + 30*traceYuAdjYu*AbsSqr(TLambdax) + 10*traceYvAdjYv*
      AbsSqr(TLambdax) + 120*AbsSqr(Lambdax)*AbsSqr(TLambdax) + 30*
      traceAdjYuTYu*Conj(TLambdax)*Lambdax + 10*traceAdjYvTYv*Conj(TLambdax)*
      Lambdax - 165*AbsSqr(MassWB)*Quad(g2) - 480*AbsSqr(MassU)*Quad(gp)*Quad(
      QHd) + 8*traceYdAdjYd*AbsSqr(MassB)*Sqr(g1) - 24*traceYeAdjYe*AbsSqr(
      MassB)*Sqr(g1) + 12*traceAdjYeTYe*Conj(MassB)*Sqr(g1) - 18*AbsSqr(MassWB)
      *Sqr(g1)*Sqr(g2) - 9*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) - 320*
      traceYdAdjYd*AbsSqr(MassG)*Sqr(g3) + 160*traceAdjYdTYd*Conj(MassG)*Sqr(g3
      ) + 72*Qd*QHd*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) + 72*Qe*QHd*AbsSqr(MassU)*Sqr
      (g1)*Sqr(gp) + 24*QHd*QHu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 72*QHd*Ql*
      AbsSqr(MassU)*Sqr(g1)*Sqr(gp) + 72*QHd*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) -
      144*QHd*Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) + 36*MassB*Qd*QHd*Conj(MassU)*
      Sqr(g1)*Sqr(gp) + 36*MassB*Qe*QHd*Conj(MassU)*Sqr(g1)*Sqr(gp) + 12*MassB*
      QHd*QHu*Conj(MassU)*Sqr(g1)*Sqr(gp) - 36*MassB*QHd*Ql*Conj(MassU)*Sqr(g1)
      *Sqr(gp) + 36*MassB*QHd*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp) - 72*MassB*QHd*Qu*
      Conj(MassU)*Sqr(g1)*Sqr(gp) - 120*traceYdAdjYd*AbsSqr(MassU)*Sqr(gp)*Sqr(
      Qd) + 60*traceAdjYdTYd*Conj(MassU)*Sqr(gp)*Sqr(Qd) - 40*traceYeAdjYe*
      AbsSqr(MassU)*Sqr(gp)*Sqr(Qe) + 20*traceAdjYeTYe*Conj(MassU)*Sqr(gp)*Sqr(
      Qe) + 120*traceYdAdjYd*AbsSqr(MassU)*Sqr(gp)*Sqr(QHd) + 40*traceYeAdjYe*
      AbsSqr(MassU)*Sqr(gp)*Sqr(QHd) + 20*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd)
      + 20*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 20*ms2*AbsSqr(Lambdax)*Sqr(
      gp)*Sqr(QHd) + 40*AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHd) + 20*
      AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHd) - 60*traceAdjYdTYd*Conj(MassU)*Sqr(gp)*
      Sqr(QHd) - 20*traceAdjYeTYe*Conj(MassU)*Sqr(gp)*Sqr(QHd) - 20*MassU*Conj(
      TLambdax)*Lambdax*Sqr(gp)*Sqr(QHd) - 48*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*Sqr
      (QHd) - 24*MassB*Conj(MassU)*Sqr(g1)*Sqr(gp)*Sqr(QHd) - 120*AbsSqr(MassU)
      *Sqr(g2)*Sqr(gp)*Sqr(QHd) - 120*AbsSqr(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(QHd) -
      60*MassWB*Conj(MassU)*Sqr(g2)*Sqr(gp)*Sqr(QHd) - 60*MassU*Conj(MassWB)*
      Sqr(g2)*Sqr(gp)*Sqr(QHd) - 1080*AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr(QHd) -
      360*AbsSqr(MassU)*Quad(gp)*Sqr(Qe)*Sqr(QHd) - 20*mHd2*AbsSqr(Lambdax)*Sqr
      (gp)*Sqr(QHu) - 20*mHu2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(QHu) - 20*ms2*AbsSqr(
      Lambdax)*Sqr(gp)*Sqr(QHu) - 40*AbsSqr(MassU)*AbsSqr(Lambdax)*Sqr(gp)*Sqr(
      QHu) - 20*AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHu) + 20*MassU*Conj(TLambdax)*
      Lambdax*Sqr(gp)*Sqr(QHu) - 240*AbsSqr(MassU)*Quad(gp)*Sqr(QHd)*Sqr(QHu) -
      40*traceYeAdjYe*AbsSqr(MassU)*Sqr(gp)*Sqr(Ql) + 20*traceAdjYeTYe*Conj(
      MassU)*Sqr(gp)*Sqr(Ql) - 720*AbsSqr(MassU)*Quad(gp)*Sqr(QHd)*Sqr(Ql) -
      120*traceYdAdjYd*AbsSqr(MassU)*Sqr(gp)*Sqr(Qq) + 60*traceAdjYdTYd*Conj(
      MassU)*Sqr(gp)*Sqr(Qq) - 2160*AbsSqr(MassU)*Quad(gp)*Sqr(QHd)*Sqr(Qq) -
      20*mHd2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) - 20*mHu2*AbsSqr(Lambdax)*Sqr(gp)
      *Sqr(Qs) - 20*ms2*AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) - 40*AbsSqr(MassU)*
      AbsSqr(Lambdax)*Sqr(gp)*Sqr(Qs) - 20*AbsSqr(TLambdax)*Sqr(gp)*Sqr(Qs) +
      20*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(Qs) - 120*AbsSqr(MassU)*Quad(
      gp)*Sqr(QHd)*Sqr(Qs) - 1080*AbsSqr(MassU)*Quad(gp)*Sqr(QHd)*Sqr(Qu) - 360
      *AbsSqr(MassU)*Quad(gp)*Sqr(QHd)*Sqr(Qv) + 60*mHd2*Sqr(Conj(Lambdax))*Sqr
      (Lambdax) + 60*mHu2*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 60*ms2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) + 30*traceconjTYuTpYu*Conj(Lambdax)*TLambdax + 10*
      traceconjTYvTpYv*Conj(Lambdax)*TLambdax - 20*Conj(MassU)*Conj(Lambdax)*
      Sqr(gp)*Sqr(QHd)*TLambdax + 20*Conj(MassU)*Conj(Lambdax)*Sqr(gp)*Sqr(QHu)
      *TLambdax + 20*Conj(MassU)*Conj(Lambdax)*Sqr(gp)*Sqr(Qs)*TLambdax));

   beta_mHd2 = beta_mHd2_1 + beta_mHd2_2;


   return beta_mHd2;
}

/**
 * Calculates the 3-loop beta function of mHd2.
 *
 * @return 3-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

/**
 * Calculates the 4-loop beta function of mHd2.
 *
 * @return 4-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

/**
 * Calculates the 5-loop beta function of mHd2.
 *
 * @return 5-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_mHd2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

} // namespace flexiblesusy
