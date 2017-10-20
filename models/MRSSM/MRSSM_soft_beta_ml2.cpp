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
 * Calculates the 1-loop beta function of ml2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_ml2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (oneOver16PiSqr*(2*mHd2*(Ye.adjoint()*Ye) + ml2*Ye.adjoint(
      )*Ye + 2*(Ye.adjoint()*me2*Ye) + Ye.adjoint()*Ye*ml2 - 0.7745966692414834
      *g1*Tr11*UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the 2-loop beta function of ml2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_ml2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (twoLoop*(0.2*(-10*(2*mHd2 + mRd2 + mS2)*AbsSqr(LamSD) - 15
      *(2*mHd2 + mRd2 + mT2)*AbsSqr(LamTD) + 2*(-5*(3*tracemd2YdAdjYd +
      traceme2YeAdjYe + traceml2AdjYeYe + 3*tracemq2AdjYdYd + 6*mHd2*
      traceYdAdjYd + 2*mHd2*traceYeAdjYe) + 6*mHd2*Sqr(g1)))*(Ye.adjoint()*Ye)
      + (-3*traceYdAdjYd - traceYeAdjYe - AbsSqr(LamSD) - 1.5*AbsSqr(LamTD) +
      1.2*Sqr(g1))*(ml2*Ye.adjoint()*Ye) + (-6*traceYdAdjYd - 2*traceYeAdjYe -
      2*AbsSqr(LamSD) - 3*AbsSqr(LamTD) + 2.4*Sqr(g1))*(Ye.adjoint()*me2*Ye) +
      (-3*traceYdAdjYd - traceYeAdjYe - AbsSqr(LamSD) - 1.5*AbsSqr(LamTD) + 1.2
      *Sqr(g1))*(Ye.adjoint()*Ye*ml2) - 8*mHd2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye
      ) - 2*(ml2*Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 4*(Ye.adjoint()*me2*Ye*
      Ye.adjoint()*Ye) - 4*(Ye.adjoint()*Ye*ml2*Ye.adjoint()*Ye) - 4*(
      Ye.adjoint()*Ye*Ye.adjoint()*me2*Ye) - 2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye
      *ml2) + (0.4*g1*(3*g1*Tr2U111 - 7.745966692414834*Tr31) + 6*Tr22*Quad(g2)
      )*UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the 3-loop beta function of ml2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> MRSSM_soft_parameters::calc_beta_ml2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

} // namespace flexiblesusy
