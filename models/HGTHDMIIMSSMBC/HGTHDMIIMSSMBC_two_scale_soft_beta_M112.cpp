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

// File generated at Tue 12 Jul 2016 10:37:46

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
 * Calculates the one-loop beta function of M112.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M112_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_M112;

   beta_M112 = Re(oneOver16PiSqr*(12*Lambda1*M112 - 12*Lambda6*M122 + 4*
      Lambda3*M222 + 2*Lambda4*M222 + 6*M112*traceYdAdjYd + 2*M112*traceYeAdjYe
      - 0.9*M112*Sqr(g1) + 3*M112*Sqr(g1d) - 6*AbsSqr(MassWB)*Sqr(g1d) + M112*
      Sqr(g1dp) - 2*AbsSqr(MassB)*Sqr(g1dp) - 4.5*M112*Sqr(g2) - 6*Sqr(g1d)*Sqr
      (Mu) - 2*Sqr(g1dp)*Sqr(Mu)));


   return beta_M112;
}

/**
 * Calculates the two-loop beta function of M112.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M112_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_M112;

   beta_M112 = Re(twoLoop*(4.6425*Power(g1,4)*M112 - 5.625*Power(g1d,4)*
      M112 - 1.125*Power(g1dp,4)*M112 - 0.1875*Power(g2,4)*M112 - 2*Lambda3*
      Lambda4*M112 + 72*Lambda1*Lambda6*M122 + 24*Lambda3*Lambda6*M122 + 24*
      Lambda4*Lambda6*M122 + 24*Lambda5*Lambda6*M122 + 12*Lambda3*Lambda7*M122
      + 12*Lambda4*Lambda7*M122 + 12*Lambda5*Lambda7*M122 + 0.9*Power(g1,4)*
      M222 + 7.5*Power(g2,4)*M222 - 8*Lambda3*Lambda4*M222 - 72*Lambda1*M112*
      traceYdAdjYd + 36*Lambda6*M122*traceYdAdjYd - 13.5*M112*
      traceYdAdjYdYdAdjYd - 4.5*M112*traceYdAdjYuYuAdjYd - 24*Lambda1*M112*
      traceYeAdjYe + 12*Lambda6*M122*traceYeAdjYe - 4.5*M112*
      traceYeAdjYeYeAdjYe + 36*Lambda6*M122*traceYuAdjYu - 24*Lambda3*M222*
      traceYuAdjYu - 12*Lambda4*M222*traceYuAdjYu + 6*g1dp*g2up*Lambda6*MassB*
      Mu + 18*g1d*g2u*Lambda6*MassWB*Mu + 14.4*Lambda1*M112*Sqr(g1) - 14.4*
      Lambda6*M122*Sqr(g1) + 4.8*Lambda3*M222*Sqr(g1) + 2.4*Lambda4*M222*Sqr(g1
      ) + 1.25*M112*traceYdAdjYd*Sqr(g1) + 3.75*M112*traceYeAdjYe*Sqr(g1) - 36*
      Lambda1*M112*Sqr(g1d) + 18*Lambda6*M122*Sqr(g1d) + 1.125*M112*Sqr(g1)*Sqr
      (g1d) - 12*Lambda1*M112*Sqr(g1dp) + 6*Lambda6*M122*Sqr(g1dp) + 0.375*M112
      *Sqr(g1)*Sqr(g1dp) - 2.25*M112*Sqr(g1d)*Sqr(g1dp) + 72*Lambda1*M112*Sqr(
      g2) - 72*Lambda6*M122*Sqr(g2) + 24*Lambda3*M222*Sqr(g2) + 12*Lambda4*M222
      *Sqr(g2) + 11.25*M112*traceYdAdjYd*Sqr(g2) + 3.75*M112*traceYeAdjYe*Sqr(
      g2) + 1.125*M112*Sqr(g1)*Sqr(g2) + 20.625*M112*Sqr(g1d)*Sqr(g2) + 1.875*
      M112*Sqr(g1dp)*Sqr(g2) + 18*Lambda6*M122*Sqr(g2u) - 12*Lambda3*M222*Sqr(
      g2u) - 6*Lambda4*M222*Sqr(g2u) - 2.25*M112*Sqr(g1d)*Sqr(g2u) + 1.5*Conj(
      MassWB)*(13*Power(g1d,4)*MassWB - 24*Power(g2,4)*MassWB + 12*g1d*g2u*
      Lambda6*Mu + Sqr(g1d)*((2*MassB + 3*MassWB)*Sqr(g1dp) + 6*MassWB*(-4*Sqr(
      g2) + Sqr(g2u)))) + 6*Lambda6*M122*Sqr(g2up) - 4*Lambda3*M222*Sqr(g2up) -
      2*Lambda4*M222*Sqr(g2up) - 0.75*M112*Sqr(g1dp)*Sqr(g2up) + 0.5*g1dp*Conj
      (MassB)*(11*Power(g1dp,3)*MassB + 12*g2up*Lambda6*Mu + g1dp*(9*MassB + 6*
      MassWB)*Sqr(g1d) + 6*g1dp*MassB*Sqr(g2up)) + 40*M112*traceYdAdjYd*Sqr(g3)
      - 60*M112*Sqr(Lambda1) - 2*M112*Sqr(Lambda3) - 8*M222*Sqr(Lambda3) - 2*
      M112*Sqr(Lambda4) - 8*M222*Sqr(Lambda4) - 3*M112*Sqr(Lambda5) - 12*M222*
      Sqr(Lambda5) - 27*M112*Sqr(Lambda6) - 18*M222*Sqr(Lambda6) + 3*M112*Sqr(
      Lambda7) - 18*M222*Sqr(Lambda7) - 2.16*Power(g1,4)*Sqr(Mu) + 18*Power(g1d
      ,4)*Sqr(Mu) + 4*Power(g1dp,4)*Sqr(Mu) - 18*Power(g2,4)*Sqr(Mu) + 6*g1d*
      g1dp*g2u*g2up*Sqr(Mu) - 1.8*Sqr(g1)*Sqr(g1d)*Sqr(Mu) - 0.6*Sqr(g1)*Sqr(
      g1dp)*Sqr(Mu) + 6*Sqr(g1d)*Sqr(g1dp)*Sqr(Mu) - 21*Sqr(g1d)*Sqr(g2)*Sqr(Mu
      ) - 3*Sqr(g1dp)*Sqr(g2)*Sqr(Mu) + 10.5*Sqr(g1d)*Sqr(g2u)*Sqr(Mu) + 1.5*
      Sqr(g1dp)*Sqr(g2u)*Sqr(Mu) + 1.5*Sqr(g1d)*Sqr(g2up)*Sqr(Mu) + 4.5*Sqr(
      g1dp)*Sqr(g2up)*Sqr(Mu)));


   return beta_M112;
}

/**
 * Calculates the three-loop beta function of M112.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M112_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

} // namespace flexiblesusy
