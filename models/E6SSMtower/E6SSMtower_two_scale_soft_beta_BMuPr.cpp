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

// File generated at Mon 19 Sep 2016 09:40:34

#include "E6SSMtower_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of BMuPr.
 *
 * @return one-loop beta function
 */
double E6SSMtower_soft_parameters::calc_beta_BMuPr_one_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPr;

   beta_BMuPr = Re(0.2*oneOver16PiSqr*(-(BMuPr*(3*Sqr(g1) + 15*Sqr(g2) +
      2*Sqr(gN))) + 2*MuPr*(3*MassB*Sqr(g1) + 15*MassWB*Sqr(g2) + 2*MassBp*Sqr(
      gN))));


   return beta_BMuPr;
}

/**
 * Calculates the two-loop beta function of BMuPr.
 *
 * @return two-loop beta function
 */
double E6SSMtower_soft_parameters::calc_beta_BMuPr_two_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPr;

   beta_BMuPr = Re(-0.06*twoLoop*(-(BMuPr*(99*Power(g1,4) + 275*Power(g2,
      4) + 64*Power(gN,4) + 20*Sqr(g2)*Sqr(gN) + 6*Sqr(g1)*(5*Sqr(g2) + 2*Sqr(
      gN)))) + 4*MuPr*(99*Power(g1,4)*MassB + 64*Power(gN,4)*MassBp + 275*Power
      (g2,4)*MassWB + 10*(MassBp + MassWB)*Sqr(g2)*Sqr(gN) + 3*Sqr(g1)*(5*(
      MassB + MassWB)*Sqr(g2) + 2*(MassB + MassBp)*Sqr(gN)))));


   return beta_BMuPr;
}

/**
 * Calculates the three-loop beta function of BMuPr.
 *
 * @return three-loop beta function
 */
double E6SSMtower_soft_parameters::calc_beta_BMuPr_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuPr;

   beta_BMuPr = 0;


   return beta_BMuPr;
}

} // namespace flexiblesusy