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

// File generated at Fri 8 Jan 2016 12:30:17

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
 * Calculates the one-loop beta function of ms2.
 *
 * @return one-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_ms2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qs = INPUT(Qs);
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_ms2;

   beta_ms2 = Re(oneOver16PiSqr*(2*gp*Qs*Tr14 + 4*(mHd2 + mHu2 + ms2)*
      AbsSqr(Lambdax) + 4*AbsSqr(TLambdax) - 8*AbsSqr(MassU)*Sqr(gp)*Sqr(Qs)));


   return beta_ms2;
}

/**
 * Calculates the two-loop beta function of ms2.
 *
 * @return two-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_ms2_two_loop(const Soft_traces& soft_traces) const
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
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_ms2;

   beta_ms2 = Re(twoLoop*(8*gp*Qs*Tr34 - 12*traceYdAdjYd*AbsSqr(TLambdax)
      - 4*traceYeAdjYe*AbsSqr(TLambdax) - 12*traceYuAdjYu*AbsSqr(TLambdax) - 4
      *traceYvAdjYv*AbsSqr(TLambdax) - 12*traceAdjYdTYd*Conj(TLambdax)*Lambdax
      - 4*traceAdjYeTYe*Conj(TLambdax)*Lambdax - 12*traceAdjYuTYu*Conj(TLambdax
      )*Lambdax - 4*traceAdjYvTYv*Conj(TLambdax)*Lambdax + 2.4*AbsSqr(TLambdax)
      *Sqr(g1) - 2.4*MassB*Conj(TLambdax)*Lambdax*Sqr(g1) + 12*AbsSqr(TLambdax)
      *Sqr(g2) - 12*MassWB*Conj(TLambdax)*Lambdax*Sqr(g2) + 8*AbsSqr(TLambdax)*
      Sqr(gp)*Sqr(QHd) - 8*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(QHd) + 8*
      AbsSqr(TLambdax)*Sqr(gp)*Sqr(QHu) - 8*MassU*Conj(TLambdax)*Lambdax*Sqr(gp
      )*Sqr(QHu) + 8*Tr2U144*Sqr(gp)*Sqr(Qs) - 8*AbsSqr(TLambdax)*Sqr(gp)*Sqr(
      Qs) + 8*MassU*Conj(TLambdax)*Lambdax*Sqr(gp)*Sqr(Qs) - 16*(mHd2 + mHu2 +
      ms2)*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 8*Conj(MassU)*Sqr(gp)*(3*MassU*Sqr
      (gp)*Sqr(Qs)*(9*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql)
      + 18*Sqr(Qq) + 3*Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv)) + Conj(Lambdax)*(Sqr(
      QHd) + Sqr(QHu) - Sqr(Qs))*(2*MassU*Lambdax - TLambdax)) + 0.8*Conj(
      Lambdax)*(-15*traceconjTYdTpTYd*Lambdax - 5*traceconjTYeTpTYe*Lambdax -
      15*traceconjTYuTpTYu*Lambdax - 5*traceconjTYvTpTYv*Lambdax - 15*
      tracemd2YdAdjYd*Lambdax - 5*traceme2YeAdjYe*Lambdax - 5*traceml2AdjYeYe*
      Lambdax - 15*tracemq2AdjYdYd*Lambdax - 15*tracemq2AdjYuYu*Lambdax - 15*
      tracemu2YuAdjYu*Lambdax - 30*mHd2*traceYdAdjYd*Lambdax - 15*mHu2*
      traceYdAdjYd*Lambdax - 15*ms2*traceYdAdjYd*Lambdax - 10*mHd2*traceYeAdjYe
      *Lambdax - 5*mHu2*traceYeAdjYe*Lambdax - 5*ms2*traceYeAdjYe*Lambdax - 15*
      mHd2*traceYuAdjYu*Lambdax - 30*mHu2*traceYuAdjYu*Lambdax - 15*ms2*
      traceYuAdjYu*Lambdax - 5*mHd2*traceYvAdjYv*Lambdax - 10*mHu2*traceYvAdjYv
      *Lambdax - 5*ms2*traceYvAdjYv*Lambdax - 5*traceYvAdjYvconjml2*Lambdax - 5
      *traceYvconjmvR2AdjYv*Lambdax - 40*AbsSqr(TLambdax)*Lambdax + 3*mHd2*
      Lambdax*Sqr(g1) + 3*mHu2*Lambdax*Sqr(g1) + 3*ms2*Lambdax*Sqr(g1) + 15*
      mHd2*Lambdax*Sqr(g2) + 15*mHu2*Lambdax*Sqr(g2) + 15*ms2*Lambdax*Sqr(g2) +
      10*mHd2*Lambdax*Sqr(gp)*Sqr(QHd) + 10*mHu2*Lambdax*Sqr(gp)*Sqr(QHd) + 10
      *ms2*Lambdax*Sqr(gp)*Sqr(QHd) + 10*mHd2*Lambdax*Sqr(gp)*Sqr(QHu) + 10*
      mHu2*Lambdax*Sqr(gp)*Sqr(QHu) + 10*ms2*Lambdax*Sqr(gp)*Sqr(QHu) - 10*mHd2
      *Lambdax*Sqr(gp)*Sqr(Qs) - 10*mHu2*Lambdax*Sqr(gp)*Sqr(Qs) - 10*ms2*
      Lambdax*Sqr(gp)*Sqr(Qs) + 3*Conj(MassB)*Sqr(g1)*(2*MassB*Lambdax -
      TLambdax) + 15*Conj(MassWB)*Sqr(g2)*(2*MassWB*Lambdax - TLambdax) - 15*
      traceconjTYdTpYd*TLambdax - 5*traceconjTYeTpYe*TLambdax - 15*
      traceconjTYuTpYu*TLambdax - 5*traceconjTYvTpYv*TLambdax)));


   return beta_ms2;
}

/**
 * Calculates the three-loop beta function of ms2.
 *
 * @return three-loop beta function
 */
double UMSSM_soft_parameters::calc_beta_ms2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_ms2;

   beta_ms2 = 0;


   return beta_ms2;
}

} // namespace flexiblesusy
