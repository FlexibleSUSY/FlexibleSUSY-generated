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

// File generated at Fri 20 Oct 2017 08:35:39

#include "HGTHDMIIMSSMBC_soft_parameters.hpp"
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
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v1_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_v1;

   beta_v1 = Re(0.1*oneOver16PiSqr*v1*(6*Sqr(g1) - 5*(6*traceYdAdjYd + 2*
      traceYeAdjYe + 3*Sqr(g1d) + Sqr(g1dp) - 6*Sqr(g2))));


   return beta_v1;
}

/**
 * Calculates the 2-loop beta function of v1.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v1_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_v1;

   beta_v1 = Re(0.00125*twoLoop*(-1407*v1*Quad(g1) - 10*v1*Sqr(g1)*(122*
      traceYdAdjYd + 174*traceYeAdjYe + 81*Sqr(g1d) + 27*Sqr(g1dp) - 45*Sqr(g2)
      ) + 25*(-32*Lambda3*Lambda4*v1 + 216*traceYdAdjYdYdAdjYd*v1 + 72*
      traceYdAdjYuYuAdjYd*v1 + 72*traceYeAdjYeYeAdjYe*v1 - 96*Lambda1*Lambda6*
      v2 - 48*Lambda3*Lambda6*v2 - 48*Lambda4*Lambda6*v2 - 48*Lambda5*Lambda6*
      v2 - 96*Lambda2*Lambda7*v2 - 48*Lambda3*Lambda7*v2 - 48*Lambda4*Lambda7*
      v2 - 48*Lambda5*Lambda7*v2 + 90*v1*Quad(g1d) + 18*v1*Quad(g1dp) + 237*v1*
      Quad(g2) - 324*traceYdAdjYd*v1*Sqr(g2) - 108*traceYeAdjYe*v1*Sqr(g2) + 6*
      v1*Sqr(g1d)*(6*Sqr(g1dp) - 67*Sqr(g2) + 6*Sqr(g2u)) - 6*v1*Sqr(g1dp)*(9*
      Sqr(g2) - 2*Sqr(g2up)) - 640*traceYdAdjYd*v1*Sqr(g3) - 192*v1*Sqr(Lambda1
      ) - 32*v1*Sqr(Lambda3) - 32*v1*Sqr(Lambda4) - 48*v1*Sqr(Lambda5) - 144*v1
      *Sqr(Lambda6) - 48*v1*Sqr(Lambda7))));


   return beta_v1;
}

/**
 * Calculates the 3-loop beta function of v1.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_v1_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v1;

   beta_v1 = 0;


   return beta_v1;
}

} // namespace flexiblesusy
