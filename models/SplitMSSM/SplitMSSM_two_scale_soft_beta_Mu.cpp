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

// File generated at Thu 15 Dec 2016 12:42:04

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
 * Calculates the one-loop beta function of Mu.
 *
 * @return one-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_Mu_one_loop(const Soft_traces& soft_traces) const
{


   double beta_Mu;

   beta_Mu = Re(0.05*oneOver16PiSqr*(20*gYd*gYu*Conj(MassB) + 60*g2d*g2u*
      Conj(MassWB) + Mu*(-18*Sqr(g1) + 5*(-18*Sqr(g2) + 3*Sqr(g2d) + 3*Sqr(g2u)
      + Sqr(gYd) + Sqr(gYu)))));


   return beta_Mu;
}

/**
 * Calculates the two-loop beta function of Mu.
 *
 * @return two-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_Mu_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_Mu;

   beta_Mu = Re(0.00125*twoLoop*(240*g2d*g2u*Conj(MassWB)*(9*Sqr(g1) + 5*
      (29*Sqr(g2) - 2*(Sqr(g2d) + Sqr(g2u)))) + 80*gYd*gYu*Conj(MassB)*(9*Sqr(
      g1) + 5*(9*Sqr(g2) - 2*(Sqr(gYd) + Sqr(gYu)))) + Mu*(2718*Power(g1,4) +
      15*Sqr(g1)*(-36*Sqr(g2) + 11*(3*Sqr(g2d) + 3*Sqr(g2u) + Sqr(gYd) + Sqr(
      gYu))) - 25*(842*Power(g2,4) - 33*Sqr(g2)*(11*Sqr(g2d) + 11*Sqr(g2u) +
      Sqr(gYd) + Sqr(gYu)) + 4*(15*Power(g2d,4) + 15*Power(g2u,4) + 2*Power(gYd
      ,4) - 24*g2d*g2u*gYd*gYu + 2*Power(gYu,4) + 9*traceYdAdjYd*Sqr(gYd) + 3*
      traceYeAdjYe*Sqr(gYd) + 9*traceYuAdjYu*Sqr(gYd) + 9*traceYdAdjYd*Sqr(gYu)
      + 3*traceYeAdjYe*Sqr(gYu) + 9*traceYuAdjYu*Sqr(gYu) + 16*Sqr(gYd)*Sqr(
      gYu) + 9*Sqr(g2u)*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu + Sqr(
      gYd) + Sqr(gYu)) + 9*Sqr(g2d)*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu + 10*Sqr(g2u) + Sqr(gYd) + Sqr(gYu)))))));


   return beta_Mu;
}

/**
 * Calculates the three-loop beta function of Mu.
 *
 * @return three-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_Mu_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_Mu;

   beta_Mu = 0;


   return beta_Mu;
}

} // namespace flexiblesusy
