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

// File generated at Thu 15 Dec 2016 12:40:48

#include "HTHDMIIMSSMBC_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of M122.
 *
 * @return one-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_M122_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M122;

   beta_M122 = Re(oneOver16PiSqr*(-6*Lambda6*M112 + 2*Lambda3*M122 + 4*
      Lambda4*M122 + 6*Lambda5*M122 - 6*Lambda7*M222 + 3*M122*traceYdAdjYd +
      M122*traceYeAdjYe + 3*M122*traceYuAdjYu - 0.9*M122*Sqr(g1) - 4.5*M122*Sqr
      (g2)));


   return beta_M122;
}

/**
 * Calculates the two-loop beta function of M122.
 *
 * @return two-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_M122_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M122;

   beta_M122 = Re(twoLoop*(36*Lambda1*Lambda6*M112 + 12*Lambda3*Lambda6*
      M112 + 12*Lambda4*Lambda6*M112 + 12*Lambda5*Lambda6*M112 + 6*Lambda3*
      Lambda7*M112 + 6*Lambda4*Lambda7*M112 + 6*Lambda5*Lambda7*M112 + 3.7425*
      Power(g1,4)*M122 - 12.6875*Power(g2,4)*M122 - 12*Lambda1*Lambda3*M122 -
      12*Lambda2*Lambda3*M122 - 12*Lambda1*Lambda4*M122 - 12*Lambda2*Lambda4*
      M122 - 6*Lambda3*Lambda4*M122 - 12*Lambda1*Lambda5*M122 - 12*Lambda2*
      Lambda5*M122 - 12*Lambda3*Lambda5*M122 - 12*Lambda4*Lambda5*M122 - 36*
      Lambda6*Lambda7*M122 + 6*Lambda3*Lambda6*M222 + 6*Lambda4*Lambda6*M222 +
      6*Lambda5*Lambda6*M222 + 36*Lambda2*Lambda7*M222 + 12*Lambda3*Lambda7*
      M222 + 12*Lambda4*Lambda7*M222 + 12*Lambda5*Lambda7*M222 + 36*Lambda6*
      M112*traceYdAdjYd - 6*Lambda3*M122*traceYdAdjYd - 12*Lambda4*M122*
      traceYdAdjYd - 18*Lambda5*M122*traceYdAdjYd - 6.75*M122*
      traceYdAdjYdYdAdjYd - 16.5*M122*traceYdAdjYuYuAdjYd + 12*Lambda6*M112*
      traceYeAdjYe - 2*Lambda3*M122*traceYeAdjYe - 4*Lambda4*M122*traceYeAdjYe
      - 6*Lambda5*M122*traceYeAdjYe - 2.25*M122*traceYeAdjYeYeAdjYe - 6*Lambda3
      *M122*traceYuAdjYu - 12*Lambda4*M122*traceYuAdjYu - 18*Lambda5*M122*
      traceYuAdjYu + 36*Lambda7*M222*traceYuAdjYu - 6.75*M122*
      traceYuAdjYuYuAdjYu + 0.375*(-96*Lambda6*M112 + 32*Lambda3*M122 + 64*
      Lambda4*M122 + 96*Lambda5*M122 - 96*Lambda7*M222 + 15*M122*traceYdAdjYd +
      5*M122*traceYeAdjYe + 15*M122*traceYuAdjYu)*Sqr(g2) - 0.025*Sqr(g1)*(288
      *Lambda6*M112 - 96*Lambda3*M122 - 192*Lambda4*M122 - 288*Lambda5*M122 +
      288*Lambda7*M222 - 25*M122*traceYdAdjYd - 75*M122*traceYeAdjYe - 85*M122*
      traceYuAdjYu - 45*M122*Sqr(g2)) + 20*M122*traceYdAdjYd*Sqr(g3) + 20*M122*
      traceYuAdjYu*Sqr(g3) + 6*M122*Sqr(Lambda1) + 6*M122*Sqr(Lambda2) + 3*M122
      *Sqr(Lambda5) - 12*M122*Sqr(Lambda6) - 12*M122*Sqr(Lambda7)));


   return beta_M122;
}

/**
 * Calculates the three-loop beta function of M122.
 *
 * @return three-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_M122_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M122;

   beta_M122 = 0;


   return beta_M122;
}

} // namespace flexiblesusy
