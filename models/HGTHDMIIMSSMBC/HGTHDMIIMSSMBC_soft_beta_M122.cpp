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

// File generated at Sun 26 Aug 2018 14:06:45

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
 * Calculates the 1-loop beta function of M122.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M122_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M122;

   beta_M122 = Re(oneOver16PiSqr*(-6*Lambda6*M112 + 2*Lambda3*M122 + 4*Lambda4*
      M122 - 6*Lambda7*M222 + 3*M122*traceYdAdjYd + M122*traceYeAdjYe + 3*M122*
      traceYuAdjYu + 6*Conj(Lambda5)*Conj(M122) + 2*g1dp*g2up*MassB*Mu + 6*g1d*
      g2u*MassWB*Mu - 0.9*M122*Sqr(g1) + 1.5*M122*Sqr(g1d) + 0.5*M122*Sqr(g1dp)
      - 4.5*M122*Sqr(g2) + 1.5*M122*Sqr(g2u) + 0.5*M122*Sqr(g2up)));


   return beta_M122;
}

/**
 * Calculates the 2-loop beta function of M122.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M122_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M122;

   const double beta_M122_1 = Re(twoLoop*(36*Lambda1*Lambda6*M112 + 12*Lambda3*
      Lambda6*M112 + 12*Lambda4*Lambda6*M112 + 6*Lambda3*Lambda7*M112 + 6*
      Lambda4*Lambda7*M112 + 6*g1d*g1dp*g2u*g2up*M122 - 12*Lambda1*Lambda3*M122
       - 12*Lambda2*Lambda3*M122 - 12*Lambda1*Lambda4*M122 - 12*Lambda2*Lambda4
      *M122 - 6*Lambda3*Lambda4*M122 + 6*Lambda3*Lambda6*M222 + 6*Lambda4*
      Lambda6*M222 + 36*Lambda2*Lambda7*M222 + 12*Lambda3*Lambda7*M222 + 12*
      Lambda4*Lambda7*M222 + 36*Lambda6*M112*traceYdAdjYd - 6*Lambda3*M122*
      traceYdAdjYd - 12*Lambda4*M122*traceYdAdjYd - 6.75*M122*
      traceYdAdjYdYdAdjYd - 16.5*M122*traceYdAdjYuYuAdjYd + 12*Lambda6*M112*
      traceYeAdjYe - 2*Lambda3*M122*traceYeAdjYe - 4*Lambda4*M122*traceYeAdjYe
      - 2.25*M122*traceYeAdjYeYeAdjYe - 6*Lambda3*M122*traceYuAdjYu - 12*
      Lambda4*M122*traceYuAdjYu + 36*Lambda7*M222*traceYuAdjYu - 6.75*M122*
      traceYuAdjYuYuAdjYu + 3*M122*AbsSqr(Lambda5) - 12*Lambda7*M122*Conj(
      Lambda6) - 2*g1dp*g2up*Lambda3*MassB*Mu - 4*g1dp*g2up*Lambda4*MassB*Mu -
      6*g1d*g2u*Lambda3*MassWB*Mu - 12*g1d*g2u*Lambda4*MassWB*Mu - 7.5*g2u*
      MassWB*Cube(g1d)*Mu - 3.5*g2up*MassB*Cube(g1dp)*Mu - 7.5*g1d*MassWB*Cube(
      g2u)*Mu - 3.5*g1dp*MassB*Cube(g2up)*Mu + 3.7425*M122*Quad(g1) - 2.8125*
      M122*Quad(g1d) - 0.5625*M122*Quad(g1dp) - 7.6875*M122*Quad(g2) - 2.8125*
      M122*Quad(g2u) - 0.5625*M122*Quad(g2up) - 7.2*Lambda6*M112*Sqr(g1) + 2.4*
      Lambda3*M122*Sqr(g1) + 4.8*Lambda4*M122*Sqr(g1) - 7.2*Lambda7*M222*Sqr(g1
      ) + 0.625*M122*traceYdAdjYd*Sqr(g1) + 1.875*M122*traceYeAdjYe*Sqr(g1) +
      2.125*M122*traceYuAdjYu*Sqr(g1) + 0.3*g1dp*g2up*MassB*Mu*Sqr(g1) + 0.9*
      g1d*g2u*MassWB*Mu*Sqr(g1) + 18*Lambda6*M112*Sqr(g1d) - 3*Lambda3*M122*Sqr
      (g1d) - 6*Lambda4*M122*Sqr(g1d) - 1.5*g1dp*g2up*MassB*Mu*Sqr(g1d) - 3*
      g1dp*g2up*MassWB*Mu*Sqr(g1d) + 0.5625*M122*Sqr(g1)*Sqr(g1d) + 6*Lambda6*
      M112*Sqr(g1dp) - Lambda3*M122*Sqr(g1dp) - 2*Lambda4*M122*Sqr(g1dp) - 3*
      g1d*g2u*MassB*Mu*Sqr(g1dp) - 1.5*g1d*g2u*MassWB*Mu*Sqr(g1dp) + 0.1875*
      M122*Sqr(g1)*Sqr(g1dp) - 1.125*M122*Sqr(g1d)*Sqr(g1dp) - 36*Lambda6*M112*
      Sqr(g2) + 12*Lambda3*M122*Sqr(g2) + 24*Lambda4*M122*Sqr(g2) - 36*Lambda7*
      M222*Sqr(g2) + 5.625*M122*traceYdAdjYd*Sqr(g2) + 1.875*M122*traceYeAdjYe*
      Sqr(g2) + 5.625*M122*traceYuAdjYu*Sqr(g2) + 1.5*g1dp*g2up*MassB*Mu*Sqr(g2
      ) + 28.5*g1d*g2u*MassWB*Mu*Sqr(g2) + 1.125*M122*Sqr(g1)*Sqr(g2) + 10.3125
      *M122*Sqr(g1d)*Sqr(g2) + 0.9375*M122*Sqr(g1dp)*Sqr(g2) - 3*Lambda3*M122*
      Sqr(g2u) - 6*Lambda4*M122*Sqr(g2u) + 18*Lambda7*M222*Sqr(g2u) - 1.5*g1dp*
      g2up*MassB*Mu*Sqr(g2u) - 3*g1dp*g2up*MassWB*Mu*Sqr(g2u) + 0.5625*M122*Sqr
      (g1)*Sqr(g2u) - 5.25*M122*Sqr(g1d)*Sqr(g2u) + 10.3125*M122*Sqr(g2)*Sqr(
      g2u) - Lambda3*M122*Sqr(g2up) - 2*Lambda4*M122*Sqr(g2up) + 6*Lambda7*M222
      *Sqr(g2up) - 3*g1d*g2u*MassB*Mu*Sqr(g2up) - 1.5*g1d*g2u*MassWB*Mu*Sqr(
      g2up) + 0.1875*M122*Sqr(g1)*Sqr(g2up) + 0.25*M122*Sqr(g1dp)*Sqr(g2up) +
      0.9375*M122*Sqr(g2)*Sqr(g2up) - 1.125*M122*Sqr(g2u)*Sqr(g2up) + 20*M122*
      traceYdAdjYd*Sqr(g3) + 20*M122*traceYuAdjYu*Sqr(g3) + 6*M122*Sqr(Lambda1)
      + 6*M122*Sqr(Lambda2)));
   const double beta_M122_2 = Re(-0.6*twoLoop*(Conj(Lambda5)*(-10*(2*M112 +
      M222)*Conj(Lambda6) - 10*(M112 + 2*M222)*Conj(Lambda7) + 20*Lambda1*Conj(
      M122) + 20*Lambda2*Conj(M122) + 20*Lambda3*Conj(M122) + 20*Lambda4*Conj(
      M122) + 30*traceYdAdjYd*Conj(M122) + 10*traceYeAdjYe*Conj(M122) + 30*
      traceYuAdjYu*Conj(M122) + 10*g1dp*g2up*Conj(MassB)*Mu + 30*g1d*g2u*Conj(
      MassWB)*Mu - 12*Conj(M122)*Sqr(g1) + 15*Conj(M122)*Sqr(g1d) + 5*Conj(M122
      )*Sqr(g1dp) - 60*Conj(M122)*Sqr(g2) + 15*Conj(M122)*Sqr(g2u) + 5*Conj(
      M122)*Sqr(g2up)) + 20*(Lambda6*M122*Conj(Lambda7) + Conj(M122)*(Lambda6*
      Lambda7 + Sqr(Lambda6) + Sqr(Lambda7)))));

   beta_M122 = beta_M122_1 + beta_M122_2;


   return beta_M122;
}

/**
 * Calculates the 3-loop beta function of M122.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M122_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M122;

   beta_M122 = 0;


   return beta_M122;
}

/**
 * Calculates the 4-loop beta function of M122.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M122_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M122;

   beta_M122 = 0;


   return beta_M122;
}

} // namespace flexiblesusy
