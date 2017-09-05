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

// File generated at Tue 5 Sep 2017 10:32:29

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
 * Calculates the one-loop beta function of MassWB.
 *
 * @return one-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassWB;

   beta_MassWB = Re(oneOver16PiSqr*(4*g2d*g2u*Mu - 12*MassWB*Sqr(g2) +
      MassWB*Sqr(g2d) + MassWB*Sqr(g2u)));


   return beta_MassWB;
}

/**
 * Calculates the two-loop beta function of MassWB.
 *
 * @return two-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MassWB;

   beta_MassWB = Re(0.008333333333333333*twoLoop*(-9320*Power(g2,4)*
      MassWB - 435*Power(g2d,4)*MassWB - 435*Power(g2u,4)*MassWB - 360*Power(
      g2d,3)*g2u*Mu + 3*Sqr(g2u)*(MassWB*(51*Sqr(g1) + 5*(91*Sqr(g2) - 6*(6*
      traceYdAdjYd + 2*traceYeAdjYe + 6*traceYuAdjYu + Sqr(gYd)))) + 5*(8*MassB
      - 7*MassWB)*Sqr(gYu)) - 24*g2d*g2u*(10*gYd*gYu*MassB + 5*Mu*Sqr(gYd) +
      Mu*(-24*Sqr(g1) + 5*(-48*Sqr(g2) + 3*Sqr(g2u) + Sqr(gYu)))) + 3*Sqr(g2d)*
      (5*(8*MassB - 7*MassWB)*Sqr(gYd) + MassWB*(51*Sqr(g1) + 5*(91*Sqr(g2) - 6
      *(6*traceYdAdjYd + 2*traceYeAdjYe + 6*traceYuAdjYu + 14*Sqr(g2u) + Sqr(
      gYu)))))));


   return beta_MassWB;
}

/**
 * Calculates the three-loop beta function of MassWB.
 *
 * @return three-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return beta_MassWB;
}

} // namespace flexiblesusy
