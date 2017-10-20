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

// File generated at Fri 20 Oct 2017 08:45:29

#include "MRSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mS2.
 *
 * @return 1-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mS2_1_loop(const Soft_traces& soft_traces) const
{


   double beta_mS2;

   beta_mS2 = Re(4*oneOver16PiSqr*((mHd2 + mRd2 + mS2)*AbsSqr(LamSD) + (
      mHu2 + mRu2 + mS2)*AbsSqr(LamSU)));


   return beta_mS2;
}

/**
 * Calculates the 2-loop beta function of mS2.
 *
 * @return 2-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mS2_2_loop(const Soft_traces& soft_traces) const
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

   beta_mS2 = Re(-0.8*twoLoop*(AbsSqr(LamSD)*(15*(2*mHd2 + 2*mRd2 + mS2 +
      mT2)*AbsSqr(LamTD) - 3*(mHd2 + mRd2 + mS2)*Sqr(g1) + 5*(3*
      tracemd2YdAdjYd + traceme2YeAdjYe + traceml2AdjYeYe + 3*tracemq2AdjYdYd +
      6*mHd2*traceYdAdjYd + 3*mRd2*traceYdAdjYd + 3*mS2*traceYdAdjYd + 2*mHd2*
      traceYeAdjYe + mRd2*traceYeAdjYe + mS2*traceYeAdjYe - 3*(mHd2 + mRd2 +
      mS2)*Sqr(g2))) + AbsSqr(LamSU)*(20*(mHu2 + mRu2 + mS2)*AbsSqr(LamSU) - 3*
      (-5*(2*mHu2 + 2*mRu2 + mS2 + mT2)*AbsSqr(LamTU) + (mHu2 + mRu2 + mS2)*Sqr
      (g1) + 5*(-tracemq2AdjYuYu - tracemu2YuAdjYu - 2*mHu2*traceYuAdjYu - mRu2
      *traceYuAdjYu - mS2*traceYuAdjYu + (mHu2 + mRu2 + mS2)*Sqr(g2)))) + 20*(
      mHd2 + mRd2 + mS2)*Sqr(LamSD)*Sqr(Conj(LamSD))));


   return beta_mS2;
}

/**
 * Calculates the 3-loop beta function of mS2.
 *
 * @return 3-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mS2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mS2;

   beta_mS2 = 0;


   return beta_mS2;
}

} // namespace flexiblesusy
