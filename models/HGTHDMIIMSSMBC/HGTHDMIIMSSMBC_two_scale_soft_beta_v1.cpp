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

// File generated at Mon 27 Feb 2017 13:22:50

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

   beta_v1 = Re(0.1*oneOver16PiSqr*v1*(6*Sqr(g1) - 5*(6*traceYdAdjYd + 2*
      traceYeAdjYe + 3*Sqr(g1d) + Sqr(g1dp) - 6*Sqr(g2))));


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

   beta_v1 = Re(0.00125*twoLoop*(-1443*Power(g1,4)*v1 - 10*v1*Sqr(g1)*(86
      *traceYdAdjYd + 162*traceYeAdjYe + 63*Sqr(g1d) + 21*Sqr(g1dp) - 9*Sqr(g2)
      ) + 25*(90*Power(g1d,4)*v1 + 18*Power(g1dp,4)*v1 + 201*Power(g2,4)*v1 -
      32*Lambda3*Lambda4*v1 + 216*traceYdAdjYdYdAdjYd*v1 + 72*
      traceYdAdjYuYuAdjYd*v1 + 72*traceYeAdjYeYeAdjYe*v1 - 96*Lambda1*Lambda6*
      v2 - 48*Lambda3*Lambda6*v2 - 48*Lambda4*Lambda6*v2 - 48*Lambda5*Lambda6*
      v2 - 96*Lambda2*Lambda7*v2 - 48*Lambda3*Lambda7*v2 - 48*Lambda4*Lambda7*
      v2 - 48*Lambda5*Lambda7*v2 - 252*traceYdAdjYd*v1*Sqr(g2) - 84*
      traceYeAdjYe*v1*Sqr(g2) + 6*v1*Sqr(g1d)*(6*Sqr(g1dp) - 61*Sqr(g2) + 6*Sqr
      (g2u)) - 6*v1*Sqr(g1dp)*(7*Sqr(g2) - 2*Sqr(g2up)) - 640*traceYdAdjYd*v1*
      Sqr(g3) - 192*v1*Sqr(Lambda1) - 32*v1*Sqr(Lambda3) - 32*v1*Sqr(Lambda4) -
      48*v1*Sqr(Lambda5) - 144*v1*Sqr(Lambda6) - 48*v1*Sqr(Lambda7))));


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
