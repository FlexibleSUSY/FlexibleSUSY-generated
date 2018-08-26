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

// File generated at Sun 26 Aug 2018 13:47:00

#include "CE6SSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of BMuPr.
 *
 * @return 1-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_BMuPr_1_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPr;

   beta_BMuPr = Re(0.2*oneOver16PiSqr*(-(BMuPr*(3*Sqr(g1) + 15*Sqr(g2) + 2*Sqr(
      gN))) + 2*MuPr*(3*MassB*Sqr(g1) + 15*MassWB*Sqr(g2) + 2*MassBp*Sqr(gN))))
      ;


   return beta_BMuPr;
}

/**
 * Calculates the 2-loop beta function of BMuPr.
 *
 * @return 2-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_BMuPr_2_loop(const Soft_traces& soft_traces) const
{


   double beta_BMuPr;

   beta_BMuPr = Re(-0.06*twoLoop*(-(BMuPr*(99*Quad(g1) + 275*Quad(g2) + 64*Quad
      (gN) + 20*Sqr(g2)*Sqr(gN) + 6*Sqr(g1)*(5*Sqr(g2) + 2*Sqr(gN)))) + 4*MuPr*
      (99*MassB*Quad(g1) + 275*MassWB*Quad(g2) + 64*MassBp*Quad(gN) + 10*(
      MassBp + MassWB)*Sqr(g2)*Sqr(gN) + 3*Sqr(g1)*(5*(MassB + MassWB)*Sqr(g2)
      + 2*(MassB + MassBp)*Sqr(gN)))));


   return beta_BMuPr;
}

/**
 * Calculates the 3-loop beta function of BMuPr.
 *
 * @return 3-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_BMuPr_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuPr;

   beta_BMuPr = 0;


   return beta_BMuPr;
}

/**
 * Calculates the 4-loop beta function of BMuPr.
 *
 * @return 4-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_BMuPr_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_BMuPr;

   beta_BMuPr = 0;


   return beta_BMuPr;
}

} // namespace flexiblesusy
