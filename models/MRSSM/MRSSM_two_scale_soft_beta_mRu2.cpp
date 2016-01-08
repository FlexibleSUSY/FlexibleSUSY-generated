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

// File generated at Fri 8 Jan 2016 15:09:53

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
 * Calculates the one-loop beta function of mRu2.
 *
 * @return one-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRu2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mRu2;

   beta_mRu2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 + 2*(mHu2 +
      mRu2 + mS2)*AbsSqr(LamSU) + 3*(mHu2 + mRu2 + mT2)*AbsSqr(LamTU)));


   return beta_mRu2;
}

/**
 * Calculates the two-loop beta function of mRu2.
 *
 * @return two-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRu2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mRu2;

   beta_mRu2 = Re(twoLoop*(6*Power(g2,4)*Tr22 - 3.0983866769659336*g1*
      Tr31 - 4*(mHd2 + mHu2 + mRd2 + mRu2 + 2*mS2)*AbsSqr(LamSD)*AbsSqr(LamSU)
      - 9*tracemq2AdjYuYu*AbsSqr(LamTU) - 9*tracemu2YuAdjYu*AbsSqr(LamTU) - 18*
      mHu2*traceYuAdjYu*AbsSqr(LamTU) - 9*mRu2*traceYuAdjYu*AbsSqr(LamTU) - 9*
      mT2*traceYuAdjYu*AbsSqr(LamTU) - 3*mHd2*AbsSqr(LamTD)*AbsSqr(LamTU) - 3*
      mHu2*AbsSqr(LamTD)*AbsSqr(LamTU) - 3*mRd2*AbsSqr(LamTD)*AbsSqr(LamTU) - 3
      *mRu2*AbsSqr(LamTD)*AbsSqr(LamTU) - 6*mT2*AbsSqr(LamTD)*AbsSqr(LamTU) - 6
      *AbsSqr(LamSU)*(tracemq2AdjYuYu + tracemu2YuAdjYu + (2*mHu2 + mRu2 + mS2)
      *traceYuAdjYu + (2*mHu2 + 2*mRu2 + mS2 + mT2)*AbsSqr(LamTU)) + 1.2*
      Tr2U111*Sqr(g1) + 12*mHu2*AbsSqr(LamTU)*Sqr(g2) + 12*mRu2*AbsSqr(LamTU)*
      Sqr(g2) + 12*mT2*AbsSqr(LamTU)*Sqr(g2) - 12*(mHu2 + mRu2 + mS2)*Sqr(LamSU
      )*Sqr(Conj(LamSU)) - 15*mHu2*Sqr(LamTU)*Sqr(Conj(LamTU)) - 15*mRu2*Sqr(
      LamTU)*Sqr(Conj(LamTU)) - 15*mT2*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_mRu2;
}

/**
 * Calculates the three-loop beta function of mRu2.
 *
 * @return three-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mRu2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRu2;

   beta_mRu2 = 0;


   return beta_mRu2;
}

} // namespace flexiblesusy
