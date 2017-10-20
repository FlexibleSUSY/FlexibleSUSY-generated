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

// File generated at Fri 20 Oct 2017 08:45:24

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
 * Calculates the 1-loop beta function of mq2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_mq2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = (oneOver16PiSqr*(2*mHd2*(Yd.adjoint()*Yd) + 2*mHu2*(
      Yu.adjoint()*Yu) + mq2*Yd.adjoint()*Yd + mq2*Yu.adjoint()*Yu + 2*(
      Yd.adjoint()*md2*Yd) + Yd.adjoint()*Yd*mq2 + 2*(Yu.adjoint()*mu2*Yu) +
      Yu.adjoint()*Yu*mq2 + 0.2581988897471611*g1*Tr11*UNITMATRIX(3))).real();


   return beta_mq2;
}

/**
 * Calculates the 2-loop beta function of mq2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_mq2_2_loop(const Soft_traces& soft_traces) const
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
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr23 = TRACE_STRUCT.Tr23;


   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = (twoLoop*((-6*tracemd2YdAdjYd - 2*traceme2YeAdjYe - 2*
      traceml2AdjYeYe - 6*tracemq2AdjYdYd - 12*mHd2*traceYdAdjYd - 4*mHd2*
      traceYeAdjYe - 2*(2*mHd2 + mRd2 + mS2)*AbsSqr(LamSD) - 3*(2*mHd2 + mRd2 +
      mT2)*AbsSqr(LamTD) + 0.8*mHd2*Sqr(g1))*(Yd.adjoint()*Yd) + (-6*
      tracemq2AdjYuYu - 6*tracemu2YuAdjYu - 12*mHu2*traceYuAdjYu - 2*(2*mHu2 +
      mRu2 + mS2)*AbsSqr(LamSU) - 3*(2*mHu2 + mRu2 + mT2)*AbsSqr(LamTU) + 1.6*
      mHu2*Sqr(g1))*(Yu.adjoint()*Yu) + (-3*traceYdAdjYd - traceYeAdjYe -
      AbsSqr(LamSD) - 1.5*AbsSqr(LamTD) + 0.4*Sqr(g1))*(mq2*Yd.adjoint()*Yd) +
      (-3*traceYuAdjYu - AbsSqr(LamSU) - 1.5*AbsSqr(LamTU) + 0.8*Sqr(g1))*(mq2*
      Yu.adjoint()*Yu) + (-6*traceYdAdjYd - 2*traceYeAdjYe - 2*AbsSqr(LamSD) -
      3*AbsSqr(LamTD) + 0.8*Sqr(g1))*(Yd.adjoint()*md2*Yd) + (-3*traceYdAdjYd -
      traceYeAdjYe - AbsSqr(LamSD) - 1.5*AbsSqr(LamTD) + 0.4*Sqr(g1))*(
      Yd.adjoint()*Yd*mq2) + (-6*traceYuAdjYu - 2*AbsSqr(LamSU) - 3*AbsSqr(
      LamTU) + 1.6*Sqr(g1))*(Yu.adjoint()*mu2*Yu) + (-3*traceYuAdjYu - AbsSqr(
      LamSU) - 1.5*AbsSqr(LamTU) + 0.8*Sqr(g1))*(Yu.adjoint()*Yu*mq2) - 8*mHd2*
      (Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 8*mHu2*(Yu.adjoint()*Yu*Yu.adjoint()*
      Yu) - 2*(mq2*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*(mq2*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu) - 4*(Yd.adjoint()*md2*Yd*Yd.adjoint()*Yd) - 4*(
      Yd.adjoint()*Yd*mq2*Yd.adjoint()*Yd) - 4*(Yd.adjoint()*Yd*Yd.adjoint()*
      md2*Yd) - 2*(Yd.adjoint()*Yd*Yd.adjoint()*Yd*mq2) - 4*(Yu.adjoint()*mu2*
      Yu*Yu.adjoint()*Yu) - 4*(Yu.adjoint()*Yu*mq2*Yu.adjoint()*Yu) - 4*(
      Yu.adjoint()*Yu*Yu.adjoint()*mu2*Yu) - 2*(Yu.adjoint()*Yu*Yu.adjoint()*Yu
      *mq2) + 0.13333333333333333*(g1*(g1*Tr2U111 + 7.745966692414834*Tr31) +
      45*Tr22*Quad(g2) + 80*Tr23*Quad(g3))*UNITMATRIX(3))).real();


   return beta_mq2;
}

/**
 * Calculates the 3-loop beta function of mq2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_mq2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = ZEROMATRIX(3,3);


   return beta_mq2;
}

} // namespace flexiblesusy
