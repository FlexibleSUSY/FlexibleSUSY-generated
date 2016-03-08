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

// File generated at Tue 8 Mar 2016 18:05:21

#include "HTHDMIIMSSMBC_two_scale_soft_parameters.hpp"
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
double HTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_one_loop(const Soft_traces& soft_traces) const
{


   double beta_Mu;

   beta_Mu = Re(-0.9*oneOver16PiSqr*Mu*(Sqr(g1) + 5*Sqr(g2)));


   return beta_Mu;
}

/**
 * Calculates the two-loop beta function of Mu.
 *
 * @return two-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_two_loop(const Soft_traces& soft_traces) const
{


   double beta_Mu;

   beta_Mu = Re(0.0125*twoLoop*Mu*(285*Power(g1,4) - 2395*Power(g2,4) -
      54*Sqr(g1)*Sqr(g2)));


   return beta_Mu;
}

/**
 * Calculates the three-loop beta function of Mu.
 *
 * @return three-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

} // namespace flexiblesusy
