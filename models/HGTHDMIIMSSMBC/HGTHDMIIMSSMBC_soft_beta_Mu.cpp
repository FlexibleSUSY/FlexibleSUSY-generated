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

// File generated at Mon 5 Mar 2018 17:08:24

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
 * Calculates the 1-loop beta function of Mu.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_1_loop(const Soft_traces& soft_traces) const
{


   double beta_Mu;

   beta_Mu = Re(0.05*oneOver16PiSqr*Mu*(-18*Sqr(g1) + 5*(3*Sqr(g1d) + Sqr
      (g1dp) - 18*Sqr(g2) + 3*Sqr(g2u) + Sqr(g2up))));


   return beta_Mu;
}

/**
 * Calculates the 2-loop beta function of Mu.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(0.00625*twoLoop*Mu*(570*Quad(g1) + 3*Sqr(g1)*(33*Sqr(g1d)
      + 11*Sqr(g1dp) - 36*Sqr(g2) + 33*Sqr(g2u) + 11*Sqr(g2up)) - 5*(-96*g1d*
      g1dp*g2u*g2up + 60*Quad(g1d) + 8*Quad(g1dp) + 798*Quad(g2) + 60*Quad(g2u)
      + 8*Quad(g2up) + 108*traceYuAdjYu*Sqr(g2u) - 363*Sqr(g2)*Sqr(g2u) + 3*
      Sqr(g1d)*(36*traceYdAdjYd + 12*traceYeAdjYe + 12*Sqr(g1dp) - 121*Sqr(g2)
      + 20*Sqr(g2u)) + 36*traceYuAdjYu*Sqr(g2up) - 33*Sqr(g2)*Sqr(g2up) + 36*
      Sqr(g2u)*Sqr(g2up) - 3*Sqr(g1dp)*(-12*traceYdAdjYd - 4*traceYeAdjYe + 11*
      Sqr(g2) + 4*Sqr(g2up)))));


   return beta_Mu;
}

/**
 * Calculates the 3-loop beta function of Mu.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_3_loop(const Soft_traces& soft_traces) const
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
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_Mu_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

} // namespace flexiblesusy
