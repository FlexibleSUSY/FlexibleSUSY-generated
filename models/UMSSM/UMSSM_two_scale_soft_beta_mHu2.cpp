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

// File generated at Sun 28 Aug 2016 15:08:51

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
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
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

   const double beta_mHu2_1 = Re(0.04*twoLoop*(Conj(MassB)*Sqr(g1)*(621*
      MassB*Sqr(g1) + 5*(9*(2*MassB + MassWB)*Sqr(g2) + 4*(-2*traceAdjYuTYu + 4
      *MassB*traceYuAdjYu + 3*(2*MassB + MassU)*QHu*(3*Qd + 3*Qe - QHd + 2*QHu
      - 3*Ql + 3*Qq - 6*Qu)*Sqr(gp)))) + 10*(15*Power(g2,4)*Tr22 +
      7.745966692414834*g1*gp*QHu*Tr2U114 + 7.745966692414834*g1*gp*QHu*Tr2U141
      + 7.745966692414834*g1*Tr31 + 20*gp*QHu*Tr34 - 5*
      traceAdjYeTYeconjTYvTpYv - 5*traceAdjYvTpYeconjTYeTYv - 15*
      tracemd2YdAdjYuYuAdjYd - 15*tracemq2AdjYdYdAdjYuYu - 15*
      tracemq2AdjYuYuAdjYdYd - 90*tracemq2AdjYuYuAdjYuYu - 15*
      tracemu2YuAdjYdYdAdjYu - 90*tracemu2YuAdjYuYuAdjYu - 15*
      traceYdAdjTYuTYuAdjYd - 15*traceYdAdjYuTYuAdjTYd - 15*mHd2*
      traceYdAdjYuYuAdjYd - 15*mHu2*traceYdAdjYuYuAdjYd - 5*
      traceYeconjTYvTpTYvAdjYe - 15*traceYuAdjTYdTYdAdjYu - 90*
      traceYuAdjTYuTYuAdjYu - 15*traceYuAdjYdTYdAdjTYu - 90*
      traceYuAdjYuTYuAdjTYu - 90*mHu2*traceYuAdjYuYuAdjYu - 30*
      traceYvAdjTYvTYvAdjYv - 5*traceYvAdjYvconjml2TpYeconjYe - 15*
      traceYvAdjYvconjml2YvAdjYv - 5*traceYvAdjYvTpTYeconjTYe - 5*
      traceYvAdjYvTpYeconjme2conjYe - 5*mHd2*traceYvAdjYvTpYeconjYe - 5*mHu2*
      traceYvAdjYvTpYeconjYe - 5*traceYvAdjYvTpYeconjYeconjml2 - 30*
      traceYvAdjYvTYvAdjTYv - 30*mHu2*traceYvAdjYvYvAdjYv - 15*
      traceYvAdjYvYvAdjYvconjml2 - 15*traceYvAdjYvYvconjmvR2AdjYv - 5*
      traceYvconjmvR2AdjYvTpYeconjYe - 15*traceYvconjmvR2AdjYvYvAdjYv + 3*
      Tr2U111*Sqr(g1) + 4*traceconjTYuTpTYu*Sqr(g1) - 4*MassB*traceconjTYuTpYu*
      Sqr(g1) + 4*tracemq2AdjYuYu*Sqr(g1) + 4*tracemu2YuAdjYu*Sqr(g1) + 4*mHu2*
      traceYuAdjYu*Sqr(g1) + 80*traceconjTYuTpTYu*Sqr(g3) - 80*MassG*
      traceconjTYuTpYu*Sqr(g3) + 80*tracemq2AdjYuYu*Sqr(g3) + 80*
      tracemu2YuAdjYu*Sqr(g3) + 80*mHu2*traceYuAdjYu*Sqr(g3) - 80*(
      traceAdjYuTYu - 2*MassG*traceYuAdjYu)*Conj(MassG)*Sqr(g3) + 20*Tr2U144*
      Sqr(gp)*Sqr(QHu) - 30*traceconjTYuTpTYu*Sqr(gp)*Sqr(QHu) + 30*MassU*
      traceconjTYuTpYu*Sqr(gp)*Sqr(QHu) - 10*traceconjTYvTpTYv*Sqr(gp)*Sqr(QHu)
      + 10*MassU*traceconjTYvTpYv*Sqr(gp)*Sqr(QHu) - 30*tracemq2AdjYuYu*Sqr(gp
      )*Sqr(QHu) - 30*tracemu2YuAdjYu*Sqr(gp)*Sqr(QHu) - 30*mHu2*traceYuAdjYu*
      Sqr(gp)*Sqr(QHu) - 10*mHu2*traceYvAdjYv*Sqr(gp)*Sqr(QHu) - 10*
      traceYvAdjYvconjml2*Sqr(gp)*Sqr(QHu) - 10*traceYvconjmvR2AdjYv*Sqr(gp)*
      Sqr(QHu) + 10*traceconjTYvTpTYv*Sqr(gp)*Sqr(Ql) - 10*MassU*
      traceconjTYvTpYv*Sqr(gp)*Sqr(Ql) + 10*mHu2*traceYvAdjYv*Sqr(gp)*Sqr(Ql) +
      10*traceYvAdjYvconjml2*Sqr(gp)*Sqr(Ql) + 10*traceYvconjmvR2AdjYv*Sqr(gp)
      *Sqr(Ql) + 30*traceconjTYuTpTYu*Sqr(gp)*Sqr(Qq) - 30*MassU*
      traceconjTYuTpYu*Sqr(gp)*Sqr(Qq) + 30*tracemq2AdjYuYu*Sqr(gp)*Sqr(Qq) +
      30*tracemu2YuAdjYu*Sqr(gp)*Sqr(Qq) + 30*mHu2*traceYuAdjYu*Sqr(gp)*Sqr(Qq)
      + 30*traceconjTYuTpTYu*Sqr(gp)*Sqr(Qu) - 30*MassU*traceconjTYuTpYu*Sqr(
      gp)*Sqr(Qu) + 30*tracemq2AdjYuYu*Sqr(gp)*Sqr(Qu) + 30*tracemu2YuAdjYu*Sqr
      (gp)*Sqr(Qu) + 30*mHu2*traceYuAdjYu*Sqr(gp)*Sqr(Qu) + 10*
      traceconjTYvTpTYv*Sqr(gp)*Sqr(Qv) - 10*MassU*traceconjTYvTpYv*Sqr(gp)*Sqr
      (Qv) + 10*mHu2*traceYvAdjYv*Sqr(gp)*Sqr(Qv) + 10*traceYvAdjYvconjml2*Sqr(
      gp)*Sqr(Qv) + 10*traceYvconjmvR2AdjYv*Sqr(gp)*Sqr(Qv))));
   const double beta_mHu2_2 = Re(0.2*twoLoop*(3*Conj(MassWB)*Sqr(g2)*(3*(
      MassB + 2*MassWB)*Sqr(g1) + 5*(11*MassWB*Sqr(g2) + 4*(MassU + 2*MassWB)*
      Sqr(gp)*Sqr(QHu))) + 4*Conj(MassU)*Sqr(gp)*(3*(MassB + 2*MassU)*QHu*(3*Qd
      + 3*Qe - QHd + 2*QHu - 3*Ql + 3*Qq - 6*Qu)*Sqr(g1) + 5*(3*traceAdjYuTYu*
      Sqr(QHu) + traceAdjYvTYv*Sqr(QHu) - 6*MassU*traceYuAdjYu*Sqr(QHu) - 2*
      MassU*traceYvAdjYv*Sqr(QHu) + 3*(2*MassU + MassWB)*Sqr(g2)*Sqr(QHu) -
      traceAdjYvTYv*Sqr(Ql) + 2*MassU*traceYvAdjYv*Sqr(Ql) - 3*traceAdjYuTYu*
      Sqr(Qq) + 6*MassU*traceYuAdjYu*Sqr(Qq) - 3*traceAdjYuTYu*Sqr(Qu) + 6*
      MassU*traceYuAdjYu*Sqr(Qu) - traceAdjYvTYv*Sqr(Qv) + 2*MassU*traceYvAdjYv
      *Sqr(Qv) + 6*MassU*Sqr(gp)*Sqr(QHu)*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(QHd) +
      4*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv)))
      + 5*Conj(Lambdax)*(Sqr(QHd) - Sqr(QHu) + Sqr(Qs))*(2*MassU*Lambdax -
      TLambdax)) + 10*(-6*(mHd2 + mHu2 + ms2)*Sqr(Conj(Lambdax))*Sqr(Lambdax) +
      Conj(Lambdax)*(Lambdax*(-3*traceconjTYdTpTYd - traceconjTYeTpTYe - 3*
      tracemd2YdAdjYd - traceme2YeAdjYe - traceml2AdjYeYe - 3*tracemq2AdjYdYd -
      6*mHd2*traceYdAdjYd - 3*mHu2*traceYdAdjYd - 3*ms2*traceYdAdjYd - 2*mHd2*
      traceYeAdjYe - mHu2*traceYeAdjYe - ms2*traceYeAdjYe + 2*(mHd2 + mHu2 +
      ms2)*Sqr(gp)*(Sqr(QHd) - Sqr(QHu) + Sqr(Qs))) - (3*traceconjTYdTpYd +
      traceconjTYeTpYe + 12*Conj(TLambdax)*Lambdax)*TLambdax) - Conj(TLambdax)*
      (Lambdax*(3*traceAdjYdTYd + traceAdjYeTYe + 2*MassU*Sqr(gp)*(Sqr(QHd) -
      Sqr(QHu) + Sqr(Qs))) + (3*traceYdAdjYd + traceYeAdjYe - 2*Sqr(gp)*(Sqr(
      QHd) - Sqr(QHu) + Sqr(Qs)))*TLambdax))));

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
