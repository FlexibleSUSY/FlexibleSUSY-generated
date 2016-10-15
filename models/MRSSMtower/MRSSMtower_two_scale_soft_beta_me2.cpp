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

// File generated at Sat 15 Oct 2016 15:08:56

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
 * Calculates the one-loop beta function of me2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMtower_soft_parameters::calc_beta_me2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (oneOver16PiSqr*(4*mHd2*(Ye*Ye.adjoint()) + 2*(me2*Ye*
      Ye.adjoint()) + 4*(Ye*ml2*Ye.adjoint()) + 2*(Ye*Ye.adjoint()*me2) +
      1.5491933384829668*g1*Tr11*UNITMATRIX(3))).real();


   return beta_me2;
}

/**
 * Calculates the two-loop beta function of me2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMtower_soft_parameters::calc_beta_me2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (twoLoop*(-0.4*(10*(2*mHd2 + mRd2 + mS2)*AbsSqr(LamSD) + 15
      *(2*mHd2 + mRd2 + mT2)*AbsSqr(LamTD) + 2*(15*tracemd2YdAdjYd + 5*
      traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*tracemq2AdjYdYd + 30*mHd2*
      traceYdAdjYd + 10*mHd2*traceYeAdjYe + 3*mHd2*Sqr(g1) - 15*mHd2*Sqr(g2)))*
      (Ye*Ye.adjoint()) + (-6*traceYdAdjYd - 2*traceYeAdjYe - 2*AbsSqr(LamSD) -
      3*AbsSqr(LamTD) - 1.2*Sqr(g1) + 6*Sqr(g2))*(me2*Ye*Ye.adjoint()) + (-12*
      traceYdAdjYd - 4*traceYeAdjYe - 4*AbsSqr(LamSD) - 6*AbsSqr(LamTD) - 2.4*
      Sqr(g1) + 12*Sqr(g2))*(Ye*ml2*Ye.adjoint()) + (-6*traceYdAdjYd - 2*
      traceYeAdjYe - 2*AbsSqr(LamSD) - 3*AbsSqr(LamTD) - 1.2*Sqr(g1) + 6*Sqr(g2
      ))*(Ye*Ye.adjoint()*me2) - 8*mHd2*(Ye*Ye.adjoint()*Ye*Ye.adjoint()) - 2*(
      me2*Ye*Ye.adjoint()*Ye*Ye.adjoint()) - 4*(Ye*ml2*Ye.adjoint()*Ye*
      Ye.adjoint()) - 4*(Ye*Ye.adjoint()*me2*Ye*Ye.adjoint()) - 4*(Ye*
      Ye.adjoint()*Ye*ml2*Ye.adjoint()) - 2*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*
      me2) + 1.6*g1*(3*g1*Tr2U111 + 3.872983346207417*Tr31)*UNITMATRIX(3)))
      .real();


   return beta_me2;
}

/**
 * Calculates the three-loop beta function of me2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMtower_soft_parameters::calc_beta_me2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return beta_me2;
}

} // namespace flexiblesusy
