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

// File generated at Sun 26 Aug 2018 14:06:48

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
 * Calculates the 1-loop beta function of M112.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M112_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;


   double beta_M112;

   beta_M112 = Re(oneOver16PiSqr*(12*Lambda1*M112 + 4*Lambda3*M222 + 2*Lambda4*
      M222 + 6*M112*traceYdAdjYd + 2*M112*traceYeAdjYe - 6*M122*Conj(Lambda6) -
      6*Lambda6*Conj(M122) - 0.9*M112*Sqr(g1) + 3*M112*Sqr(g1d) - 6*AbsSqr(
      MassWB)*Sqr(g1d) + M112*Sqr(g1dp) - 2*AbsSqr(MassB)*Sqr(g1dp) - 4.5*M112*
      Sqr(g2) - 6*Sqr(g1d)*Sqr(Mu) - 2*Sqr(g1dp)*Sqr(Mu)));


   return beta_M112;
}

/**
 * Calculates the 2-loop beta function of M112.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M112_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;


   double beta_M112;

   const double beta_M112_1 = Re(twoLoop*(-2*Lambda3*Lambda4*M112 + 12*Lambda5*
      Lambda6*M122 + 6*Lambda5*Lambda7*M122 - 8*Lambda3*Lambda4*M222 - 72*
      Lambda1*M112*traceYdAdjYd - 13.5*M112*traceYdAdjYdYdAdjYd - 4.5*M112*
      traceYdAdjYuYuAdjYd - 24*Lambda1*M112*traceYeAdjYe - 4.5*M112*
      traceYeAdjYeYeAdjYe - 24*Lambda3*M222*traceYuAdjYu - 12*Lambda4*M222*
      traceYuAdjYu + 3*M112*AbsSqr(Lambda7) - 18*M222*AbsSqr(Lambda7) + 6*
      Lambda3*M122*Conj(Lambda7) + 6*Lambda4*M122*Conj(Lambda7) + 36*Lambda1*
      Lambda6*Conj(M122) + 12*Lambda3*Lambda6*Conj(M122) + 12*Lambda4*Lambda6*
      Conj(M122) + 6*Lambda3*Lambda7*Conj(M122) + 6*Lambda4*Lambda7*Conj(M122)
      + 18*Lambda6*traceYdAdjYd*Conj(M122) + 6*Lambda6*traceYeAdjYe*Conj(M122)
      + 18*Lambda6*traceYuAdjYu*Conj(M122) - 3*Conj(Lambda5)*(Lambda5*(M112 + 4
      *M222) - 4*Conj(Lambda6)*Conj(M122) - 2*Conj(Lambda7)*Conj(M122)) +
      4.6425*M112*Quad(g1) + 0.9*M222*Quad(g1) - 5.625*M112*Quad(g1d) - 1.125*
      M112*Quad(g1dp) - 0.1875*M112*Quad(g2) + 7.5*M222*Quad(g2) + 14.4*Lambda1
      *M112*Sqr(g1) + 4.8*Lambda3*M222*Sqr(g1) + 2.4*Lambda4*M222*Sqr(g1) +
      1.25*M112*traceYdAdjYd*Sqr(g1) + 3.75*M112*traceYeAdjYe*Sqr(g1) - 7.2*
      Lambda6*Conj(M122)*Sqr(g1) - 36*Lambda1*M112*Sqr(g1d) + 9*Lambda6*Conj(
      M122)*Sqr(g1d) + 1.125*M112*Sqr(g1)*Sqr(g1d) - 12*Lambda1*M112*Sqr(g1dp)
      + 3*Lambda6*Conj(M122)*Sqr(g1dp) + 0.375*M112*Sqr(g1)*Sqr(g1dp) - 2.25*
      M112*Sqr(g1d)*Sqr(g1dp) + 4.5*AbsSqr(MassB)*Sqr(g1d)*Sqr(g1dp) + 72*
      Lambda1*M112*Sqr(g2) + 24*Lambda3*M222*Sqr(g2) + 12*Lambda4*M222*Sqr(g2)
      + 11.25*M112*traceYdAdjYd*Sqr(g2) + 3.75*M112*traceYeAdjYe*Sqr(g2) - 36*
      Lambda6*Conj(M122)*Sqr(g2) + 1.125*M112*Sqr(g1)*Sqr(g2) + 20.625*M112*Sqr
      (g1d)*Sqr(g2) + 1.875*M112*Sqr(g1dp)*Sqr(g2) - 12*Lambda3*M222*Sqr(g2u) -
      6*Lambda4*M222*Sqr(g2u) + 9*Lambda6*Conj(M122)*Sqr(g2u) - 2.25*M112*Sqr(
      g1d)*Sqr(g2u) - 4*Lambda3*M222*Sqr(g2up) - 2*Lambda4*M222*Sqr(g2up) + 3*
      Lambda6*Conj(M122)*Sqr(g2up) - 0.75*M112*Sqr(g1dp)*Sqr(g2up) - 0.6*Conj(
      Lambda6)*(15*Lambda6*(3*M112 + 2*M222) + 12*M122*Sqr(g1) - 5*(2*g1dp*g2up
      *MassB*Mu + 6*g1d*g2u*MassWB*Mu + 3*M122*Sqr(g1d) + M122*Sqr(g1dp) + M122
      *(12*Lambda1 + 4*Lambda3 + 4*Lambda4 + 6*traceYdAdjYd + 2*traceYeAdjYe +
      6*traceYuAdjYu - 12*Sqr(g2) + 3*Sqr(g2u) + Sqr(g2up)))) + 40*M112*
      traceYdAdjYd*Sqr(g3) - 60*M112*Sqr(Lambda1) - 2*M112*Sqr(Lambda3) - 8*
      M222*Sqr(Lambda3) - 2*M112*Sqr(Lambda4) - 8*M222*Sqr(Lambda4) + 6*g1d*
      g1dp*g2u*g2up*Sqr(Mu) - 2.16*Quad(g1)*Sqr(Mu) + 18*Quad(g1d)*Sqr(Mu) + 4*
      Quad(g1dp)*Sqr(Mu) - 18*Quad(g2)*Sqr(Mu) - 1.8*Sqr(g1)*Sqr(g1d)*Sqr(Mu) -
      0.6*Sqr(g1)*Sqr(g1dp)*Sqr(Mu) + 6*Sqr(g1d)*Sqr(g1dp)*Sqr(Mu) - 21*Sqr(g1d
      )*Sqr(g2)*Sqr(Mu) - 3*Sqr(g1dp)*Sqr(g2)*Sqr(Mu) + 10.5*Sqr(g1d)*Sqr(g2u)*
      Sqr(Mu) + 1.5*Sqr(g1dp)*Sqr(g2u)*Sqr(Mu) + 1.5*Sqr(g1d)*Sqr(g2up)*Sqr(Mu)
      + 4.5*Sqr(g1dp)*Sqr(g2up)*Sqr(Mu)));
   const double beta_M112_2 = Re(0.5*twoLoop*(3*Conj(MassWB)*(12*g1d*g2u*
      Lambda6*Mu + 13*MassWB*Quad(g1d) - 24*MassWB*Quad(g2) + Sqr(g1d)*((2*
      MassB + 3*MassWB)*Sqr(g1dp) + 6*MassWB*(-4*Sqr(g2) + Sqr(g2u)))) + g1dp*
      Conj(MassB)*(11*MassB*Cube(g1dp) + 12*g2up*Lambda6*Mu + 6*g1dp*(MassWB*
      Sqr(g1d) + MassB*Sqr(g2up)))));

   beta_M112 = beta_M112_1 + beta_M112_2;


   return beta_M112;
}

/**
 * Calculates the 3-loop beta function of M112.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M112_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

/**
 * Calculates the 4-loop beta function of M112.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M112_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

} // namespace flexiblesusy
