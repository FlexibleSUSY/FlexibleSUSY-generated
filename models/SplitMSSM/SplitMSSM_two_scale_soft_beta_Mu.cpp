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

// File generated at Mon 19 Sep 2016 09:48:24

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
      (29*Sqr(g2) - 2*(Sqr(g2d) + Sqr(g2u)))) + Mu*(2718*Power(g1,4) - 21050*
      Power(g2,4) - 1500*Power(g2d,4) - 1500*Power(g2u,4) - 200*Power(gYd,4) +
      2400*g2d*g2u*gYd*gYu - 200*Power(gYu,4) - 540*Sqr(g1)*Sqr(g2) - 2700*
      traceYuAdjYu*Sqr(g2d) + 495*Sqr(g1)*Sqr(g2d) + 9075*Sqr(g2)*Sqr(g2d) -
      2700*traceYuAdjYu*Sqr(g2u) + 495*Sqr(g1)*Sqr(g2u) + 9075*Sqr(g2)*Sqr(g2u)
      - 9000*Sqr(g2d)*Sqr(g2u) - 900*traceYuAdjYu*Sqr(gYd) + 165*Sqr(g1)*Sqr(
      gYd) + 825*Sqr(g2)*Sqr(gYd) - 900*Sqr(g2d)*Sqr(gYd) - 900*Sqr(g2u)*Sqr(
      gYd) - 900*traceYuAdjYu*Sqr(gYu) + 165*Sqr(g1)*Sqr(gYu) + 825*Sqr(g2)*Sqr
      (gYu) - 900*Sqr(g2d)*Sqr(gYu) - 900*Sqr(g2u)*Sqr(gYu) - 1600*Sqr(gYd)*Sqr
      (gYu) - 900*traceYdAdjYd*(3*Sqr(g2d) + 3*Sqr(g2u) + Sqr(gYd) + Sqr(gYu))
      - 300*traceYeAdjYe*(3*Sqr(g2d) + 3*Sqr(g2u) + Sqr(gYd) + Sqr(gYu))) + 80*
      gYd*gYu*Conj(MassB)*(9*Sqr(g1) + 5*(9*Sqr(g2) - 2*(Sqr(gYd) + Sqr(gYu))))
      ));


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
