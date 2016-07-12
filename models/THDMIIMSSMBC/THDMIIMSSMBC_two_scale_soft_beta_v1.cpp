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

// File generated at Tue 12 Jul 2016 10:31:04

#include "THDMIIMSSMBC_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of v1.
 *
 * @return one-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v1_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_v1;

   beta_v1 = Re(oneOver16PiSqr*(-3*traceYdAdjYd*v1 - traceYeAdjYe*v1 +
      0.6*v1*(Sqr(g1) + 5*Sqr(g2))));


   return beta_v1;
}

/**
 * Calculates the two-loop beta function of v1.
 *
 * @return two-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v1_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_v1;

   beta_v1 = Re(twoLoop*(-1.65375*Power(g1,4)*v1 + 10.03125*Power(g2,4)*
      v1 - Lambda3*Lambda4*v1 + 6.75*traceYdAdjYdYdAdjYd*v1 + 2.25*
      traceYdAdjYuYuAdjYd*v1 + 2.25*traceYeAdjYeYeAdjYe*v1 - 3*Lambda1*Lambda6*
      v2 - 1.5*Lambda3*Lambda6*v2 - 1.5*Lambda4*Lambda6*v2 - 1.5*Lambda5*
      Lambda6*v2 - 3*Lambda2*Lambda7*v2 - 1.5*Lambda3*Lambda7*v2 - 1.5*Lambda4*
      Lambda7*v2 - 1.5*Lambda5*Lambda7*v2 + 0.1125*v1*Sqr(g1)*Sqr(g2) - 0.075*
      traceYeAdjYe*v1*(27*Sqr(g1) + 35*Sqr(g2)) - 0.025*traceYdAdjYd*v1*(43*Sqr
      (g1) + 315*Sqr(g2) + 800*Sqr(g3)) - 6*v1*Sqr(Lambda1) - v1*Sqr(Lambda3) -
      v1*Sqr(Lambda4) - 1.5*v1*Sqr(Lambda5) - 4.5*v1*Sqr(Lambda6) - 1.5*v1*Sqr
      (Lambda7)));


   return beta_v1;
}

/**
 * Calculates the three-loop beta function of v1.
 *
 * @return three-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v1_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v1;

   beta_v1 = 0;


   return beta_v1;
}

} // namespace flexiblesusy
