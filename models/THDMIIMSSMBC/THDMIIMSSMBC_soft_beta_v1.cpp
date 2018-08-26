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

// File generated at Sun 26 Aug 2018 14:08:54

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

   beta_v1 = Re(oneOver16PiSqr*(0.6*v1*Sqr(g1) + v1*(-3*traceYdAdjYd -
      traceYeAdjYe + 3*Sqr(g2))));


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

   beta_v1 = Re(twoLoop*(-(Lambda3*Lambda4*v1) + 6.75*traceYdAdjYdYdAdjYd*v1 +
      2.25*traceYdAdjYuYuAdjYd*v1 + 2.25*traceYeAdjYeYeAdjYe*v1 - 1.5*Lambda1*
      Lambda6*v2 - 0.75*Lambda3*Lambda6*v2 - 0.75*Lambda4*Lambda6*v2 - 0.75*
      Lambda5*Lambda6*v2 - 1.5*Lambda2*Lambda7*v2 - 0.75*Lambda3*Lambda7*v2 -
      0.75*Lambda4*Lambda7*v2 - 0.75*Lambda5*Lambda7*v2 - 1.5*v1*AbsSqr(Lambda7
      ) - 0.75*(6*Lambda6*v1 + (2*Lambda1 + Lambda3 + Lambda4)*v2)*Conj(Lambda6
      ) - 1.5*Lambda2*v2*Conj(Lambda7) - 0.75*Lambda3*v2*Conj(Lambda7) - 0.75*
      Lambda4*v2*Conj(Lambda7) - 0.75*Conj(Lambda5)*(2*Lambda5*v1 + v2*Conj(
      Lambda6) + v2*Conj(Lambda7)) - 1.60875*v1*Quad(g1) + 11.15625*v1*Quad(g2)
      - 1.525*traceYdAdjYd*v1*Sqr(g1) - 2.175*traceYeAdjYe*v1*Sqr(g1) - 10.125*
      traceYdAdjYd*v1*Sqr(g2) - 3.375*traceYeAdjYe*v1*Sqr(g2) + 0.5625*v1*Sqr(
      g1)*Sqr(g2) - 20*traceYdAdjYd*v1*Sqr(g3) - 6*v1*Sqr(Lambda1) - v1*Sqr(
      Lambda3) - v1*Sqr(Lambda4)));


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

} // namespace flexiblesusy
