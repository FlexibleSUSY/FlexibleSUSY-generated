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

// File generated at Sat 15 Oct 2016 15:25:20

#include "HGTHDMIIMSSMBC_two_scale_soft_parameters.hpp"
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
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v1_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_v1;

   beta_v1 = Re(0.1*oneOver16PiSqr*v1*(-30*traceYdAdjYd - 10*traceYeAdjYe
      + 6*Sqr(g1) - 15*Sqr(g1d) - 5*Sqr(g1dp) + 30*Sqr(g2)));


   return beta_v1;
}

/**
 * Calculates the two-loop beta function of v1.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v1_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_v1;

   beta_v1 = Re(0.00125*twoLoop*(-1443*Power(g1,4)*v1 + 2250*Power(g1d,4)
      *v1 + 450*Power(g1dp,4)*v1 + 5025*Power(g2,4)*v1 - 800*Lambda3*Lambda4*v1
      + 5400*traceYdAdjYdYdAdjYd*v1 + 1800*traceYdAdjYuYuAdjYd*v1 + 1800*
      traceYeAdjYeYeAdjYe*v1 - 2400*Lambda1*Lambda6*v2 - 1200*Lambda3*Lambda6*
      v2 - 1200*Lambda4*Lambda6*v2 - 1200*Lambda5*Lambda6*v2 - 2400*Lambda2*
      Lambda7*v2 - 1200*Lambda3*Lambda7*v2 - 1200*Lambda4*Lambda7*v2 - 1200*
      Lambda5*Lambda7*v2 - 630*v1*Sqr(g1)*Sqr(g1d) - 210*v1*Sqr(g1)*Sqr(g1dp) +
      900*v1*Sqr(g1d)*Sqr(g1dp) + 90*v1*Sqr(g1)*Sqr(g2) - 9150*v1*Sqr(g1d)*Sqr
      (g2) - 1050*v1*Sqr(g1dp)*Sqr(g2) - 60*traceYeAdjYe*v1*(27*Sqr(g1) + 35*
      Sqr(g2)) + 900*v1*Sqr(g1d)*Sqr(g2u) + 300*v1*Sqr(g1dp)*Sqr(g2up) - 20*
      traceYdAdjYd*v1*(43*Sqr(g1) + 315*Sqr(g2) + 800*Sqr(g3)) - 4800*v1*Sqr(
      Lambda1) - 800*v1*Sqr(Lambda3) - 800*v1*Sqr(Lambda4) - 1200*v1*Sqr(
      Lambda5) - 3600*v1*Sqr(Lambda6) - 1200*v1*Sqr(Lambda7)));


   return beta_v1;
}

/**
 * Calculates the three-loop beta function of v1.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v1_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v1;

   beta_v1 = 0;


   return beta_v1;
}

} // namespace flexiblesusy
