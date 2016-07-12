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

// File generated at Tue 12 Jul 2016 10:37:51

#include "SplitMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of v.
 *
 * @return one-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_v_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_v;

   beta_v = Re(0.1*oneOver16PiSqr*v*(-30*traceYdAdjYd - 10*traceYeAdjYe -
      30*traceYuAdjYu + 6*Sqr(g1) + 30*Sqr(g2) - 15*Sqr(g2d) - 15*Sqr(g2u) - 5
      *Sqr(gYd) - 5*Sqr(gYu)));


   return beta_v;
}

/**
 * Calculates the two-loop beta function of v.
 *
 * @return two-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_v_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_v;

   beta_v = Re(-0.00125*twoLoop*v*(1377*Power(g1,4) - 5575*Power(g2,4) -
      2250*Power(g2d,4) - 2250*Power(g2u,4) - 450*Power(gYd,4) - 2400*g2d*g2u*
      gYd*gYu - 450*Power(gYu,4) - 5400*traceYdAdjYdYdAdjYd + 1200*
      traceYdAdjYuYuAdjYd - 1800*traceYeAdjYeYeAdjYe - 5400*traceYuAdjYuYuAdjYu
      + 2060*traceYuAdjYu*Sqr(g1) + 6300*traceYuAdjYu*Sqr(g2) - 90*Sqr(g1)*Sqr
      (g2) + 60*traceYeAdjYe*(27*Sqr(g1) + 35*Sqr(g2)) + 630*Sqr(g1)*Sqr(g2d) +
      9150*Sqr(g2)*Sqr(g2d) + 630*Sqr(g1)*Sqr(g2u) + 9150*Sqr(g2)*Sqr(g2u) -
      600*Sqr(g2d)*Sqr(g2u) + 16000*traceYuAdjYu*Sqr(g3) + 20*traceYdAdjYd*(43*
      Sqr(g1) + 315*Sqr(g2) + 800*Sqr(g3)) + 210*Sqr(g1)*Sqr(gYd) + 1050*Sqr(g2
      )*Sqr(gYd) - 900*Sqr(g2d)*Sqr(gYd) + 210*Sqr(g1)*Sqr(gYu) + 1050*Sqr(g2)*
      Sqr(gYu) - 900*Sqr(g2u)*Sqr(gYu) - 1000*Sqr(gYd)*Sqr(gYu) + 1200*Sqr(
      Lambdax)));


   return beta_v;
}

/**
 * Calculates the three-loop beta function of v.
 *
 * @return three-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_v_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v;

   beta_v = 0;


   return beta_v;
}

} // namespace flexiblesusy
