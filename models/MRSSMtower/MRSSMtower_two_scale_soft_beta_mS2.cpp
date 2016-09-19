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

// File generated at Mon 19 Sep 2016 09:37:20

#include "MRSSMtower_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mS2.
 *
 * @return one-loop beta function
 */
double MRSSMtower_soft_parameters::calc_beta_mS2_one_loop(const Soft_traces& soft_traces) const
{


   double beta_mS2;

   beta_mS2 = Re(4*oneOver16PiSqr*((mHd2 + mRd2 + mS2)*AbsSqr(LamSD) + (
      mHu2 + mRu2 + mS2)*AbsSqr(LamSU)));


   return beta_mS2;
}

/**
 * Calculates the two-loop beta function of mS2.
 *
 * @return two-loop beta function
 */
double MRSSMtower_soft_parameters::calc_beta_mS2_two_loop(const Soft_traces& soft_traces) const
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


   double beta_mS2;

   beta_mS2 = Re(-0.8*twoLoop*(AbsSqr(LamSD)*(15*tracemd2YdAdjYd + 5*
      traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*tracemq2AdjYdYd + 15*(2*mHd2 +
      mRd2 + mS2)*traceYdAdjYd + 10*mHd2*traceYeAdjYe + 5*mRd2*traceYeAdjYe + 5
      *mS2*traceYeAdjYe + 15*(2*mHd2 + 2*mRd2 + mS2 + mT2)*AbsSqr(LamTD) - 3*
      mHd2*Sqr(g1) - 3*mRd2*Sqr(g1) - 3*mS2*Sqr(g1) - 15*mHd2*Sqr(g2) - 15*mRd2
      *Sqr(g2) - 15*mS2*Sqr(g2)) + AbsSqr(LamSU)*(20*(mHu2 + mRu2 + mS2)*AbsSqr
      (LamSU) - 3*(-5*tracemq2AdjYuYu - 5*tracemu2YuAdjYu - 5*(2*mHu2 + mRu2 +
      mS2)*traceYuAdjYu - 5*(2*mHu2 + 2*mRu2 + mS2 + mT2)*AbsSqr(LamTU) + mHu2*
      Sqr(g1) + mRu2*Sqr(g1) + mS2*Sqr(g1) + 5*mHu2*Sqr(g2) + 5*mRu2*Sqr(g2) +
      5*mS2*Sqr(g2))) + 20*(mHd2 + mRd2 + mS2)*Sqr(LamSD)*Sqr(Conj(LamSD))));


   return beta_mS2;
}

/**
 * Calculates the three-loop beta function of mS2.
 *
 * @return three-loop beta function
 */
double MRSSMtower_soft_parameters::calc_beta_mS2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mS2;

   beta_mS2 = 0;


   return beta_mS2;
}

} // namespace flexiblesusy
