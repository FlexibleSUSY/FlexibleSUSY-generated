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
 * Calculates the one-loop beta function of Mu.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_one_loop(const Soft_traces& soft_traces) const
{


   double beta_Mu;

   beta_Mu = Re(0.05*oneOver16PiSqr*Mu*(-18*Sqr(g1) + 5*(3*Sqr(g1d) + Sqr
      (g1dp) - 18*Sqr(g2) + 3*Sqr(g2u) + Sqr(g2up))));


   return beta_Mu;
}

/**
 * Calculates the two-loop beta function of Mu.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(0.00625*twoLoop*Mu*(570*Power(g1,4) - 300*Power(g1d,4) -
      40*Power(g1dp,4) - 3990*Power(g2,4) - 300*Power(g2u,4) + 480*g1d*g1dp*g2u
      *g2up - 40*Power(g2up,4) + 99*Sqr(g1)*Sqr(g1d) + 33*Sqr(g1)*Sqr(g1dp) -
      180*Sqr(g1d)*Sqr(g1dp) - 180*traceYdAdjYd*(3*Sqr(g1d) + Sqr(g1dp)) - 60*
      traceYeAdjYe*(3*Sqr(g1d) + Sqr(g1dp)) - 108*Sqr(g1)*Sqr(g2) + 1815*Sqr(
      g1d)*Sqr(g2) + 165*Sqr(g1dp)*Sqr(g2) - 540*traceYuAdjYu*Sqr(g2u) + 99*Sqr
      (g1)*Sqr(g2u) - 300*Sqr(g1d)*Sqr(g2u) + 1815*Sqr(g2)*Sqr(g2u) - 180*
      traceYuAdjYu*Sqr(g2up) + 33*Sqr(g1)*Sqr(g2up) + 60*Sqr(g1dp)*Sqr(g2up) +
      165*Sqr(g2)*Sqr(g2up) - 180*Sqr(g2u)*Sqr(g2up)));


   return beta_Mu;
}

/**
 * Calculates the three-loop beta function of Mu.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

} // namespace flexiblesusy
