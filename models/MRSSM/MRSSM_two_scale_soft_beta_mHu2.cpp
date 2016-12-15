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

// File generated at Thu 15 Dec 2016 12:46:16

#include "MRSSM_two_scale_soft_parameters.hpp"
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
double MRSSM_soft_parameters::calc_beta_mHu2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mHu2;

   beta_mHu2 = Re(oneOver16PiSqr*(0.7745966692414834*g1*Tr11 + 6*
      tracemq2AdjYuYu + 6*tracemu2YuAdjYu + 6*mHu2*traceYuAdjYu + 2*(mHu2 +
      mRu2 + mS2)*AbsSqr(LamSU) + 3*(mHu2 + mRu2 + mT2)*AbsSqr(LamTU)));


   return beta_mHu2;
}

/**
 * Calculates the two-loop beta function of mHu2.
 *
 * @return two-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mHu2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
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
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mHu2;

   beta_mHu2 = Re(twoLoop*(6*Power(g2,4)*Tr22 + 3.0983866769659336*g1*
      Tr31 - 6*tracemd2YdAdjYuYuAdjYd - 6*tracemq2AdjYdYdAdjYuYu - 6*
      tracemq2AdjYuYuAdjYdYd - 36*tracemq2AdjYuYuAdjYuYu - 6*
      tracemu2YuAdjYdYdAdjYu - 36*tracemu2YuAdjYuYuAdjYu - 6*mHd2*
      traceYdAdjYuYuAdjYd - 6*mHu2*traceYdAdjYuYuAdjYd - 36*mHu2*
      traceYuAdjYuYuAdjYu - 4*(mHd2 + mHu2 + mRd2 + mRu2 + 2*mS2)*AbsSqr(LamSD)
      *AbsSqr(LamSU) - 6*(2*mHu2 + 2*mRu2 + mS2 + mT2)*AbsSqr(LamSU)*AbsSqr(
      LamTU) - 3*mHd2*AbsSqr(LamTD)*AbsSqr(LamTU) - 3*mHu2*AbsSqr(LamTD)*AbsSqr
      (LamTU) - 3*mRd2*AbsSqr(LamTD)*AbsSqr(LamTU) - 3*mRu2*AbsSqr(LamTD)*
      AbsSqr(LamTU) - 6*mT2*AbsSqr(LamTD)*AbsSqr(LamTU) + 1.2*Tr2U111*Sqr(g1) +
      1.6*tracemq2AdjYuYu*Sqr(g1) + 1.6*tracemu2YuAdjYu*Sqr(g1) + 1.6*mHu2*
      traceYuAdjYu*Sqr(g1) + 12*mHu2*AbsSqr(LamTU)*Sqr(g2) + 12*mRu2*AbsSqr(
      LamTU)*Sqr(g2) + 12*mT2*AbsSqr(LamTU)*Sqr(g2) + 32*tracemq2AdjYuYu*Sqr(g3
      ) + 32*tracemu2YuAdjYu*Sqr(g3) + 32*mHu2*traceYuAdjYu*Sqr(g3) - 12*(mHu2
      + mRu2 + mS2)*Sqr(LamSU)*Sqr(Conj(LamSU)) - 15*mHu2*Sqr(LamTU)*Sqr(Conj(
      LamTU)) - 15*mRu2*Sqr(LamTU)*Sqr(Conj(LamTU)) - 15*mT2*Sqr(LamTU)*Sqr(
      Conj(LamTU))));


   return beta_mHu2;
}

/**
 * Calculates the three-loop beta function of mHu2.
 *
 * @return three-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mHu2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

} // namespace flexiblesusy
