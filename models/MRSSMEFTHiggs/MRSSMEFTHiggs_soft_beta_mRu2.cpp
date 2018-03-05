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

// File generated at Mon 5 Mar 2018 15:32:03

#include "MRSSMEFTHiggs_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mRu2.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mRu2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   double beta_mRu2;

   beta_mRu2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 + 2*(mHu2 +
      mRu2 + mS2)*AbsSqr(LamSU) + 3*(mHu2 + mRu2 + mT2)*AbsSqr(LamTU)));


   return beta_mRu2;
}

/**
 * Calculates the 2-loop beta function of mRu2.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mRu2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mRu2;

   beta_mRu2 = Re(twoLoop*(-3.0983866769659336*g1*Tr31 - 2*AbsSqr(LamSU)*
      (2*(mHd2 + mHu2 + mRd2 + mRu2 + 2*mS2)*AbsSqr(LamSD) + 3*(tracemq2AdjYuYu
      + tracemu2YuAdjYu + 2*mHu2*traceYuAdjYu + mRu2*traceYuAdjYu + mS2*
      traceYuAdjYu + (2*mHu2 + 2*mRu2 + mS2 + mT2)*AbsSqr(LamTU))) + 6*Tr22*
      Quad(g2) + 1.2*Tr2U111*Sqr(g1) + 3*AbsSqr(LamTU)*(-3*(tracemq2AdjYuYu +
      tracemu2YuAdjYu + (2*mHu2 + mRu2 + mT2)*traceYuAdjYu) - (mHd2 + mHu2 +
      mRd2 + mRu2 + 2*mT2)*AbsSqr(LamTD) + 4*(mHu2 + mRu2 + mT2)*Sqr(g2)) - 12*
      (mHu2 + mRu2 + mS2)*Sqr(LamSU)*Sqr(Conj(LamSU)) - 15*(mHu2 + mRu2 + mT2)*
      Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_mRu2;
}

/**
 * Calculates the 3-loop beta function of mRu2.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mRu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRu2;

   beta_mRu2 = 0;


   return beta_mRu2;
}

/**
 * Calculates the 4-loop beta function of mRu2.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mRu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mRu2;

   beta_mRu2 = 0;


   return beta_mRu2;
}

} // namespace flexiblesusy
