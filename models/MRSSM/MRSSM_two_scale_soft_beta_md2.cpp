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

// File generated at Sun 28 Aug 2016 15:06:30

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
 * Calculates the one-loop beta function of md2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_md2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (oneOver16PiSqr*(4*mHd2*(Yd*Yd.adjoint()) + 2*(md2*Yd*
      Yd.adjoint()) + 4*(Yd*mq2*Yd.adjoint()) + 2*(Yd*Yd.adjoint()*md2) +
      0.5163977794943222*g1*Tr11*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the two-loop beta function of md2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_md2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (twoLoop*(-0.4*(30*tracemd2YdAdjYd + 10*traceme2YeAdjYe +
      10*traceml2AdjYeYe + 30*tracemq2AdjYdYd + 60*mHd2*traceYdAdjYd + 20*mHd2*
      traceYeAdjYe + 10*(2*mHd2 + mRd2 + mS2)*AbsSqr(LamSD) + 15*(2*mHd2 + mRd2
      + mT2)*AbsSqr(LamTD) - 2*mHd2*Sqr(g1) - 30*mHd2*Sqr(g2))*(Yd*Yd.adjoint(
      )) + (-6*traceYdAdjYd - 2*traceYeAdjYe - 2*AbsSqr(LamSD) - 3*AbsSqr(LamTD
      ) + 0.4*Sqr(g1) + 6*Sqr(g2))*(md2*Yd*Yd.adjoint()) + (-12*traceYdAdjYd -
      4*traceYeAdjYe - 4*AbsSqr(LamSD) - 6*AbsSqr(LamTD) + 0.8*Sqr(g1) + 12*Sqr
      (g2))*(Yd*mq2*Yd.adjoint()) + (-6*traceYdAdjYd - 2*traceYeAdjYe - 2*
      AbsSqr(LamSD) - 3*AbsSqr(LamTD) + 0.4*Sqr(g1) + 6*Sqr(g2))*(Yd*Yd.adjoint
      ()*md2) - 8*mHd2*(Yd*Yd.adjoint()*Yd*Yd.adjoint()) + (-4*mHd2 - 4*mHu2)*(
      Yd*Yu.adjoint()*Yu*Yd.adjoint()) - 2*(md2*Yd*Yd.adjoint()*Yd*Yd.adjoint()
      ) - 2*(md2*Yd*Yu.adjoint()*Yu*Yd.adjoint()) - 4*(Yd*mq2*Yd.adjoint()*Yd*
      Yd.adjoint()) - 4*(Yd*mq2*Yu.adjoint()*Yu*Yd.adjoint()) - 4*(Yd*
      Yd.adjoint()*md2*Yd*Yd.adjoint()) - 4*(Yd*Yd.adjoint()*Yd*mq2*Yd.adjoint(
      )) - 2*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*md2) - 4*(Yd*Yu.adjoint()*mu2*Yu*
      Yd.adjoint()) - 4*(Yd*Yu.adjoint()*Yu*mq2*Yd.adjoint()) - 2*(Yd*
      Yu.adjoint()*Yu*Yd.adjoint()*md2) + 0.5333333333333333*(20*Power(g3,4)*
      Tr23 + g1*(g1*Tr2U111 + 3.872983346207417*Tr31))*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the three-loop beta function of md2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_md2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

} // namespace flexiblesusy
