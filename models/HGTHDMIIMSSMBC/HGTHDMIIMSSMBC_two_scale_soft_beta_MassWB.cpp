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

// File generated at Mon 9 May 2016 12:04:21

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
 * Calculates the one-loop beta function of MassWB.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_MassWB_one_loop(const Soft_traces& soft_traces) const
{


   double beta_MassWB;

   beta_MassWB = Re(MassWB*oneOver16PiSqr*(Sqr(g1d) - 12*Sqr(g2) + Sqr(
      g2u)));


   return beta_MassWB;
}

/**
 * Calculates the two-loop beta function of MassWB.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_MassWB_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_MassWB;

   beta_MassWB = Re(0.025*twoLoop*(-145*Power(g1d,4)*MassWB - 2960*Power(
      g2,4)*MassWB - 145*Power(g2u,4)*MassWB - 180*MassWB*traceYdAdjYd*Sqr(g1d)
      - 60*MassWB*traceYeAdjYe*Sqr(g1d) + 51*MassWB*Sqr(g1)*Sqr(g1d) + 40*
      MassB*Sqr(g1d)*Sqr(g1dp) - 35*MassWB*Sqr(g1d)*Sqr(g1dp) + 455*MassWB*Sqr(
      g1d)*Sqr(g2) - 180*MassWB*traceYuAdjYu*Sqr(g2u) + 51*MassWB*Sqr(g1)*Sqr(
      g2u) + 455*MassWB*Sqr(g2)*Sqr(g2u) + 40*MassB*Sqr(g2u)*Sqr(g2up) - 35*
      MassWB*Sqr(g2u)*Sqr(g2up)));


   return beta_MassWB;
}

/**
 * Calculates the three-loop beta function of MassWB.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_MassWB_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_MassWB;

   beta_MassWB = 0;


   return beta_MassWB;
}

} // namespace flexiblesusy
