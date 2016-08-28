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

// File generated at Sun 28 Aug 2016 15:06:32

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
 * Calculates the one-loop beta function of mT2.
 *
 * @return one-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mT2_one_loop(const Soft_traces& soft_traces) const
{


   double beta_mT2;

   beta_mT2 = Re(2*oneOver16PiSqr*((mHd2 + mRd2 + mT2)*AbsSqr(LamTD) + (
      mHu2 + mRu2 + mT2)*AbsSqr(LamTU)));


   return beta_mT2;
}

/**
 * Calculates the two-loop beta function of mT2.
 *
 * @return two-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mT2_two_loop(const Soft_traces& soft_traces) const
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

   beta_mT2 = Re(-0.4*twoLoop*(-40*Power(g2,4)*Tr22 + AbsSqr(LamTD)*(15*
      tracemd2YdAdjYd + 5*traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*
      tracemq2AdjYdYd + 15*(2*mHd2 + mRd2 + mT2)*traceYdAdjYd + 10*mHd2*
      traceYeAdjYe + 5*mRd2*traceYeAdjYe + 5*mT2*traceYeAdjYe + 10*(2*mHd2 + 2*
      mRd2 + mS2 + mT2)*AbsSqr(LamSD) - 3*mHd2*Sqr(g1) - 3*mRd2*Sqr(g1) - 3*mT2
      *Sqr(g1) + 5*mHd2*Sqr(g2) + 5*mRd2*Sqr(g2) + 5*mT2*Sqr(g2)) + AbsSqr(
      LamTU)*(15*tracemq2AdjYuYu + 15*tracemu2YuAdjYu + 15*(2*mHu2 + mRu2 + mT2
      )*traceYuAdjYu + 10*(2*mHu2 + 2*mRu2 + mS2 + mT2)*AbsSqr(LamSU) - 3*mHu2*
      Sqr(g1) - 3*mRu2*Sqr(g1) - 3*mT2*Sqr(g1) + 5*mHu2*Sqr(g2) + 5*mRu2*Sqr(g2
      ) + 5*mT2*Sqr(g2)) + 30*(mHd2 + mRd2 + mT2)*Sqr(LamTD)*Sqr(Conj(LamTD)) +
      30*(mHu2 + mRu2 + mT2)*Sqr(LamTU)*Sqr(Conj(LamTU))));


   return beta_mT2;
}

/**
 * Calculates the three-loop beta function of mT2.
 *
 * @return three-loop beta function
 */
double MRSSM_soft_parameters::calc_beta_mT2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mT2;

   beta_mT2 = 0;


   return beta_mT2;
}

} // namespace flexiblesusy
