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

// File generated at Sun 4 Aug 2019 17:36:08

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
 * Calculates the 1-loop beta function of me2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_soft_parameters::calc_beta_me2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (oneOver16PiSqr*(4*mHd2*(Ye*Ye.adjoint()) + 2*(me2*Ye*Ye.adjoint(
      )) + 4*(Ye*ml2*Ye.adjoint()) + 2*(Ye*Ye.adjoint()*me2) +
      1.5491933384829668*g1*Tr11*UNITMATRIX(3))).real();


   return beta_me2;
}

/**
 * Calculates the 2-loop beta function of me2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_soft_parameters::calc_beta_me2_2_loop(const Soft_traces& soft_traces) const
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

   beta_me2 = (twoLoop*(-0.4*(30*tracemd2YdAdjYd + 10*traceme2YeAdjYe + 10*
      traceml2AdjYeYe + 30*tracemq2AdjYdYd + 60*mHd2*traceYdAdjYd + 20*mHd2*
      traceYeAdjYe + 20*mHd2*AbsSqr(LamSD) + 10*mRd2*AbsSqr(LamSD) + 10*mS2*
      AbsSqr(LamSD) + 30*mHd2*AbsSqr(LamTD) + 15*mRd2*AbsSqr(LamTD) + 15*mT2*
      AbsSqr(LamTD) + 6*mHd2*Sqr(g1) - 30*mHd2*Sqr(g2))*(Ye*Ye.adjoint()) + 0.2
      *(-30*traceYdAdjYd - 10*traceYeAdjYe - 10*AbsSqr(LamSD) - 15*AbsSqr(LamTD
      ) - 6*Sqr(g1) + 30*Sqr(g2))*(me2*Ye*Ye.adjoint()) - 0.4*(30*traceYdAdjYd
      + 10*traceYeAdjYe + 10*AbsSqr(LamSD) + 15*AbsSqr(LamTD) + 6*Sqr(g1) - 30*
      Sqr(g2))*(Ye*ml2*Ye.adjoint()) + 0.2*(-30*traceYdAdjYd - 10*traceYeAdjYe
      - 10*AbsSqr(LamSD) - 15*AbsSqr(LamTD) - 6*Sqr(g1) + 30*Sqr(g2))*(Ye*Ye.
      adjoint()*me2) - 8*mHd2*(Ye*Ye.adjoint()*Ye*Ye.adjoint()) - 2*(me2*Ye*Ye.
      adjoint()*Ye*Ye.adjoint()) - 4*(Ye*ml2*Ye.adjoint()*Ye*Ye.adjoint()) - 4*
      (Ye*Ye.adjoint()*me2*Ye*Ye.adjoint()) - 4*(Ye*Ye.adjoint()*Ye*ml2*Ye.
      adjoint()) - 2*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*me2) + 1.6*g1*(3*g1*
      Tr2U111 + 3.872983346207417*Tr31)*UNITMATRIX(3))).real();


   return beta_me2;
}

/**
 * Calculates the 3-loop beta function of me2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_soft_parameters::calc_beta_me2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return beta_me2;
}

/**
 * Calculates the 4-loop beta function of me2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_soft_parameters::calc_beta_me2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return beta_me2;
}

/**
 * Calculates the 5-loop beta function of me2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSMEFTHiggs_soft_parameters::calc_beta_me2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return beta_me2;
}

} // namespace flexiblesusy
