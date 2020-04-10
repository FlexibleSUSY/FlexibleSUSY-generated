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

// File generated at Fri 10 Apr 2020 19:42:52

#include "THDMIIMSSMBC_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of v1.
 *
 * @return 1-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v1_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_v1;

   beta_v1 = Re(0.2*oneOver16PiSqr*v1*(-15*traceYdAdjYd - 5*traceYeAdjYe + 3*
      Sqr(g1) + 15*Sqr(g2)));


   return beta_v1;
}

/**
 * Calculates the 2-loop beta function of v1.
 *
 * @return 2-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v1_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_v1;

   beta_v1 = Re(0.00125*twoLoop*(-800*Lambda3*Lambda4*v1 + 5400*
      traceYdAdjYdYdAdjYd*v1 + 1800*traceYdAdjYuYuAdjYd*v1 + 1800*
      traceYeAdjYeYeAdjYe*v1 - 1200*Lambda1*Lambda6*v2 - 600*Lambda3*Lambda6*v2
       - 600*Lambda4*Lambda6*v2 - 600*Lambda5*Lambda6*v2 - 1200*Lambda2*Lambda7
      *v2 - 600*Lambda3*Lambda7*v2 - 600*Lambda4*Lambda7*v2 - 600*Lambda5*
      Lambda7*v2 - 1200*v1*AbsSqr(Lambda5) - 3600*v1*AbsSqr(Lambda6) - 1200*v1*
      AbsSqr(Lambda7) - 1200*Lambda1*v2*Conj(Lambda6) - 600*Lambda3*v2*Conj(
      Lambda6) - 600*Lambda4*v2*Conj(Lambda6) - 600*v2*Conj(Lambda5)*Conj(
      Lambda6) - 1200*Lambda2*v2*Conj(Lambda7) - 600*Lambda3*v2*Conj(Lambda7) -
      600*Lambda4*v2*Conj(Lambda7) - 600*v2*Conj(Lambda5)*Conj(Lambda7) - 1287*
      v1*Quad(g1) + 8925*v1*Quad(g2) - 1220*traceYdAdjYd*v1*Sqr(g1) - 1740*
      traceYeAdjYe*v1*Sqr(g1) - 8100*traceYdAdjYd*v1*Sqr(g2) - 2700*
      traceYeAdjYe*v1*Sqr(g2) + 450*v1*Sqr(g1)*Sqr(g2) - 16000*traceYdAdjYd*v1*
      Sqr(g3) - 4800*v1*Sqr(Lambda1) - 800*v1*Sqr(Lambda3) - 800*v1*Sqr(Lambda4
      )));


   return beta_v1;
}

/**
 * Calculates the 3-loop beta function of v1.
 *
 * @return 3-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v1_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v1;

   beta_v1 = 0;


   return beta_v1;
}

/**
 * Calculates the 4-loop beta function of v1.
 *
 * @return 4-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v1_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v1;

   beta_v1 = 0;


   return beta_v1;
}

/**
 * Calculates the 5-loop beta function of v1.
 *
 * @return 5-loop beta function
 */
double THDMIIMSSMBC_soft_parameters::calc_beta_v1_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v1;

   beta_v1 = 0;


   return beta_v1;
}

} // namespace flexiblesusy
