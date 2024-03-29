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


#include "MRSSM2_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mu2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_soft_parameters::calc_beta_mu2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (4*mHu2*(Yu*Yu.adjoint()) + 2*(mu2*Yu*Yu.adjoint()) + 4*(Yu*mq2*
      Yu.adjoint()) + 2*(Yu*Yu.adjoint()*mu2) - 1.0327955589886444*g1*Tr11*
      UNITMATRIX(3)).real();


   return oneLoop * beta_mu2;
}

/**
 * Calculates the 2-loop beta function of mu2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_soft_parameters::calc_beta_mu2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (-0.4*(30*tracemq2AdjYuYu + 30*tracemu2YuAdjYu + 60*mHu2*
      traceYuAdjYu + 20*mHu2*AbsSqr(LamSU) + 10*mRu2*AbsSqr(LamSU) + 10*mS2*
      AbsSqr(LamSU) + 30*mHu2*AbsSqr(LamTU) + 15*mRu2*AbsSqr(LamTU) + 15*mT2*
      AbsSqr(LamTU) + 2*mHu2*Sqr(g1) - 30*mHu2*Sqr(g2))*(Yu*Yu.adjoint()) + 0.2
      *(-30*traceYuAdjYu - 10*AbsSqr(LamSU) - 15*AbsSqr(LamTU) - 2*Sqr(g1) + 30
      *Sqr(g2))*(mu2*Yu*Yu.adjoint()) - 0.4*(30*traceYuAdjYu + 10*AbsSqr(LamSU)
      + 15*AbsSqr(LamTU) + 2*Sqr(g1) - 30*Sqr(g2))*(Yu*mq2*Yu.adjoint()) + 0.2*
      (-30*traceYuAdjYu - 10*AbsSqr(LamSU) - 15*AbsSqr(LamTU) - 2*Sqr(g1) + 30*
      Sqr(g2))*(Yu*Yu.adjoint()*mu2) - 4*(mHd2 + mHu2)*(Yu*Yd.adjoint()*Yd*Yu.
      adjoint()) - 8*mHu2*(Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 2*(mu2*Yu*Yd.
      adjoint()*Yd*Yu.adjoint()) - 2*(mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 4*
      (Yu*mq2*Yd.adjoint()*Yd*Yu.adjoint()) - 4*(Yu*mq2*Yu.adjoint()*Yu*Yu.
      adjoint()) - 4*(Yu*Yd.adjoint()*md2*Yd*Yu.adjoint()) - 4*(Yu*Yd.adjoint()
      *Yd*mq2*Yu.adjoint()) - 2*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*mu2) - 4*(Yu*
      Yu.adjoint()*mu2*Yu*Yu.adjoint()) - 4*(Yu*Yu.adjoint()*Yu*mq2*Yu.adjoint(
      )) - 2*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*mu2) + 1.0666666666666667*(-
      3.872983346207417*g1*Tr31 + 10*Tr23*Quad(g3) + 2*Tr2U111*Sqr(g1))*
      UNITMATRIX(3)).real();


   return twoLoop * beta_mu2;
}

/**
 * Calculates the 3-loop beta function of mu2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_soft_parameters::calc_beta_mu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return threeLoop * beta_mu2;
}

/**
 * Calculates the 4-loop beta function of mu2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_soft_parameters::calc_beta_mu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return fourLoop * beta_mu2;
}

/**
 * Calculates the 5-loop beta function of mu2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM2_soft_parameters::calc_beta_mu2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return fiveLoop * beta_mu2;
}

} // namespace flexiblesusy
