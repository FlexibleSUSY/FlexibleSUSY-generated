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

// File generated at Tue 22 Jan 2019 16:12:27

#include "HTHDMIIMSSMBC_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of Mu.
 *
 * @return 1-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_1_loop(const Soft_traces& soft_traces) const
{


   double beta_Mu;

   beta_Mu = Re(-0.9*oneOver16PiSqr*Mu*(Sqr(g1) + 5*Sqr(g2)));


   return beta_Mu;
}

/**
 * Calculates the 2-loop beta function of Mu.
 *
 * @return 2-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_2_loop(const Soft_traces& soft_traces) const
{


   double beta_Mu;

   beta_Mu = Re(0.0125*twoLoop*Mu*(285*Quad(g1) - 2395*Quad(g2) - 54*Sqr(g1)*
      Sqr(g2)));


   return beta_Mu;
}

/**
 * Calculates the 3-loop beta function of Mu.
 *
 * @return 3-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

/**
 * Calculates the 4-loop beta function of Mu.
 *
 * @return 4-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

/**
 * Calculates the 5-loop beta function of Mu.
 *
 * @return 5-loop beta function
 */
double HTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

} // namespace flexiblesusy
