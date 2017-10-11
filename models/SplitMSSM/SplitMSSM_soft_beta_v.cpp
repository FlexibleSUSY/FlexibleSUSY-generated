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

// File generated at Tue 10 Oct 2017 21:10:04

#include "SplitMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of v.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_v_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_v;

   beta_v = Re(0.1*oneOver16PiSqr*v*(6*Sqr(g1) + 5*(-6*traceYdAdjYd - 2*
      traceYeAdjYe - 6*traceYuAdjYu + 6*Sqr(g2) - 3*Sqr(g2d) - 3*Sqr(g2u) - Sqr
      (gYd) - Sqr(gYu))));


   return beta_v;
}

/**
 * Calculates the 2-loop beta function of v.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_v_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_v;

   beta_v = Re(-0.00125*twoLoop*v*(1341*Quad(g1) + 10*Sqr(g1)*(122*
      traceYdAdjYd + 174*traceYeAdjYe + 242*traceYuAdjYu - 45*Sqr(g2) + 81*Sqr(
      g2d) + 81*Sqr(g2u) + 27*Sqr(gYd) + 27*Sqr(gYu)) - 25*(259*Quad(g2) - 6*
      Sqr(g2)*(67*Sqr(g2d) + 67*Sqr(g2u) + 9*(6*traceYdAdjYd + 2*traceYeAdjYe +
      6*traceYuAdjYu + Sqr(gYd) + Sqr(gYu))) + 2*(48*g2d*g2u*gYd*gYu + 108*
      traceYdAdjYdYdAdjYd - 24*traceYdAdjYuYuAdjYd + 36*traceYeAdjYeYeAdjYe +
      108*traceYuAdjYuYuAdjYu + 45*Quad(g2d) + 45*Quad(g2u) + 9*Quad(gYd) + 9*
      Quad(gYu) - 320*traceYdAdjYd*Sqr(g3) - 320*traceYuAdjYu*Sqr(g3) + 6*Sqr(
      g2d)*(2*Sqr(g2u) + 3*Sqr(gYd)) + 18*Sqr(g2u)*Sqr(gYu) + 20*Sqr(gYd)*Sqr(
      gYu) - 24*Sqr(Lambdax)))));


   return beta_v;
}

/**
 * Calculates the 3-loop beta function of v.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_v_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_v;

   beta_v = 0;


   return beta_v;
}

} // namespace flexiblesusy
