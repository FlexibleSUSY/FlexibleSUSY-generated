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

// File generated at Mon 5 Mar 2018 17:42:07

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
 * Calculates the 1-loop beta function of MassWB.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_1_loop(const Soft_traces& soft_traces) const
{


   double beta_MassWB;

   beta_MassWB = Re(oneOver16PiSqr*(4*g2d*g2u*Mu - 12*MassWB*Sqr(g2) +
      MassWB*Sqr(g2d) + MassWB*Sqr(g2u)));


   return beta_MassWB;
}

/**
 * Calculates the 2-loop beta function of MassWB.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MassWB;

   beta_MassWB = Re(0.008333333333333333*twoLoop*(-360*g2u*Cube(g2d)*Mu -
      9320*MassWB*Quad(g2) - 435*MassWB*Quad(g2d) - 435*MassWB*Quad(g2u) + 3*
      Sqr(g2u)*(MassWB*(51*Sqr(g1) + 5*(91*Sqr(g2) - 6*(6*traceYdAdjYd + 2*
      traceYeAdjYe + 6*traceYuAdjYu + Sqr(gYd)))) + 5*(8*MassB - 7*MassWB)*Sqr(
      gYu)) - 24*g2d*g2u*(10*gYd*gYu*MassB + 5*Mu*Sqr(gYd) + Mu*(-24*Sqr(g1) +
      5*(-48*Sqr(g2) + 3*Sqr(g2u) + Sqr(gYu)))) + 3*Sqr(g2d)*(5*(8*MassB - 7*
      MassWB)*Sqr(gYd) + MassWB*(51*Sqr(g1) + 5*(91*Sqr(g2) - 6*(6*traceYdAdjYd
      + 2*traceYeAdjYe + 6*traceYuAdjYu + 14*Sqr(g2u) + Sqr(gYu)))))));


   return beta_MassWB;
}

/**
 * Calculates the 3-loop beta function of MassWB.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return beta_MassWB;
}

/**
 * Calculates the 4-loop beta function of MassWB.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_MassWB_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return beta_MassWB;
}

} // namespace flexiblesusy
