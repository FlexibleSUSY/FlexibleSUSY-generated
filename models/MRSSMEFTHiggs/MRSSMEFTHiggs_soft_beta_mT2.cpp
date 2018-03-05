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

// File generated at Mon 5 Mar 2018 15:32:01

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
 * Calculates the 1-loop beta function of mT2.
 *
 * @return 1-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mT2_1_loop(const Soft_traces& soft_traces) const
{


   double beta_mT2;

   beta_mT2 = Re(2*oneOver16PiSqr*((mHd2 + mRd2 + mT2)*AbsSqr(LamTD) + (
      mHu2 + mRu2 + mT2)*AbsSqr(LamTU)));


   return beta_mT2;
}

/**
 * Calculates the 2-loop beta function of mT2.
 *
 * @return 2-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mT2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr22 = TRACE_STRUCT.Tr22;


   double beta_mT2;

   beta_mT2 = Re(0.4*twoLoop*(40*Tr22*Quad(g2) + AbsSqr(LamTD)*(-10*(2*
      mHd2 + 2*mRd2 + mS2 + mT2)*AbsSqr(LamSD) + 3*(mHd2 + mRd2 + mT2)*Sqr(g1)
      - 5*(3*tracemd2YdAdjYd + traceme2YeAdjYe + traceml2AdjYeYe + 3*
      tracemq2AdjYdYd + 6*mHd2*traceYdAdjYd + 3*mRd2*traceYdAdjYd + 3*mT2*
      traceYdAdjYd + 2*mHd2*traceYeAdjYe + mRd2*traceYeAdjYe + mT2*traceYeAdjYe
      + (mHd2 + mRd2 + mT2)*Sqr(g2))) + AbsSqr(LamTU)*(-10*(2*mHu2 + 2*mRu2 +
      mS2 + mT2)*AbsSqr(LamSU) + 3*(mHu2 + mRu2 + mT2)*Sqr(g1) - 5*(3*(
      tracemq2AdjYuYu + tracemu2YuAdjYu + (2*mHu2 + mRu2 + mT2)*traceYuAdjYu) +
      (mHu2 + mRu2 + mT2)*Sqr(g2))) - 30*(mHd2 + mRd2 + mT2)*Sqr(LamTD)*Sqr(
      Conj(LamTD)) - 30*(mHu2 + mRu2 + mT2)*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_mT2;
}

/**
 * Calculates the 3-loop beta function of mT2.
 *
 * @return 3-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mT2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mT2;

   beta_mT2 = 0;


   return beta_mT2;
}

/**
 * Calculates the 4-loop beta function of mT2.
 *
 * @return 4-loop beta function
 */
double MRSSMEFTHiggs_soft_parameters::calc_beta_mT2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mT2;

   beta_mT2 = 0;


   return beta_mT2;
}

} // namespace flexiblesusy
