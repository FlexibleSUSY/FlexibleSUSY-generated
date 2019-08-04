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

// File generated at Sun 4 Aug 2019 19:38:47

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

   beta_M112 = Re(0.1*oneOver16PiSqr*(120*Lambda1*M112 + 40*Lambda3*M222 + 20*
      Lambda4*M222 + 60*M112*traceYdAdjYd + 20*M112*traceYeAdjYe - 60*M122*Conj
      (Lambda6) - 60*Lambda6*Conj(M122) - 9*M112*Sqr(g1) + 30*M112*Sqr(g1d) -
      60*AbsSqr(MassWB)*Sqr(g1d) + 10*M112*Sqr(g1dp) - 20*AbsSqr(MassB)*Sqr(
      g1dp) - 45*M112*Sqr(g2) - 60*Sqr(g1d)*Sqr(Mu) - 20*Sqr(g1dp)*Sqr(Mu)));


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

   const double beta_M112_1 = Re(0.0025*twoLoop*(-800*Lambda3*Lambda4*M112 +
      4200*Lambda5*Lambda6*M122 + 1800*Lambda5*Lambda7*M122 - 3200*Lambda3*
      Lambda4*M222 - 28800*Lambda1*M112*traceYdAdjYd - 5400*M112*
      traceYdAdjYdYdAdjYd - 1800*M112*traceYdAdjYuYuAdjYd - 9600*Lambda1*M112*
      traceYeAdjYe - 1800*M112*traceYeAdjYeYeAdjYe - 9600*Lambda3*M222*
      traceYuAdjYu - 4800*Lambda4*M222*traceYuAdjYu - 1200*M112*AbsSqr(Lambda5)
      - 4800*M222*AbsSqr(Lambda5) - 10800*M112*AbsSqr(Lambda6) - 7200*M222*
      AbsSqr(Lambda6) + 1200*M112*AbsSqr(Lambda7) - 7200*M222*AbsSqr(Lambda7) +
      13200*Lambda1*M122*Conj(Lambda6) + 4200*Lambda3*M122*Conj(Lambda6) + 4200
      *Lambda4*M122*Conj(Lambda6) + 7200*M122*traceYdAdjYd*Conj(Lambda6) + 2400
      *M122*traceYeAdjYe*Conj(Lambda6) + 7200*M122*traceYuAdjYu*Conj(Lambda6) -
      1200*Lambda2*M122*Conj(Lambda7) + 1800*Lambda3*M122*Conj(Lambda7) + 1800*
      Lambda4*M122*Conj(Lambda7) + 13200*Lambda1*Lambda6*Conj(M122) + 4200*
      Lambda3*Lambda6*Conj(M122) + 4200*Lambda4*Lambda6*Conj(M122) - 1200*
      Lambda2*Lambda7*Conj(M122) + 1800*Lambda3*Lambda7*Conj(M122) + 1800*
      Lambda4*Lambda7*Conj(M122) + 7200*Lambda6*traceYdAdjYd*Conj(M122) + 2400*
      Lambda6*traceYeAdjYe*Conj(M122) + 7200*Lambda6*traceYuAdjYu*Conj(M122) +
      4200*Conj(Lambda5)*Conj(Lambda6)*Conj(M122) + 2400*g1dp*g2up*MassB*Conj(
      Lambda6)*Mu + 7200*g1d*g2u*MassWB*Conj(Lambda6)*Mu + 1857*M112*Quad(g1) +
      360*M222*Quad(g1) - 2250*M112*Quad(g1d) - 450*M112*Quad(g1dp) - 75*M112*
      Quad(g2) + 3000*M222*Quad(g2) + 5760*Lambda1*M112*Sqr(g1) + 1920*Lambda3*
      M222*Sqr(g1) + 960*Lambda4*M222*Sqr(g1) + 500*M112*traceYdAdjYd*Sqr(g1) +
      1500*M112*traceYeAdjYe*Sqr(g1) - 2880*M122*Conj(Lambda6)*Sqr(g1) - 2880*
      Lambda6*Conj(M122)*Sqr(g1) - 14400*Lambda1*M112*Sqr(g1d) + 3600*M122*Conj
      (Lambda6)*Sqr(g1d) + 3600*Lambda6*Conj(M122)*Sqr(g1d) + 450*M112*Sqr(g1)*
      Sqr(g1d) - 4800*Lambda1*M112*Sqr(g1dp) + 1200*M122*Conj(Lambda6)*Sqr(g1dp
      ) + 1200*Lambda6*Conj(M122)*Sqr(g1dp) + 150*M112*Sqr(g1)*Sqr(g1dp) - 900*
      M112*Sqr(g1d)*Sqr(g1dp) + 28800*Lambda1*M112*Sqr(g2) + 9600*Lambda3*M222*
      Sqr(g2) + 4800*Lambda4*M222*Sqr(g2) + 4500*M112*traceYdAdjYd*Sqr(g2) +
      1500*M112*traceYeAdjYe*Sqr(g2) - 14400*M122*Conj(Lambda6)*Sqr(g2) - 14400
      *Lambda6*Conj(M122)*Sqr(g2) + 450*M112*Sqr(g1)*Sqr(g2) + 8250*M112*Sqr(
      g1d)*Sqr(g2) + 750*M112*Sqr(g1dp)*Sqr(g2) - 4800*Lambda3*M222*Sqr(g2u) -
      2400*Lambda4*M222*Sqr(g2u) + 3600*M122*Conj(Lambda6)*Sqr(g2u) + 3600*
      Lambda6*Conj(M122)*Sqr(g2u) - 900*M112*Sqr(g1d)*Sqr(g2u) - 1600*Lambda3*
      M222*Sqr(g2up) - 800*Lambda4*M222*Sqr(g2up) + 1200*M122*Conj(Lambda6)*Sqr
      (g2up) + 1200*Lambda6*Conj(M122)*Sqr(g2up) - 300*M112*Sqr(g1dp)*Sqr(g2up)
      + 16000*M112*traceYdAdjYd*Sqr(g3) - 24000*M112*Sqr(Lambda1) - 800*M112*
      Sqr(Lambda3) - 3200*M222*Sqr(Lambda3) - 800*M112*Sqr(Lambda4) - 3200*M222
      *Sqr(Lambda4) + 2400*g1d*g1dp*g2u*g2up*Sqr(Mu) - 864*Quad(g1)*Sqr(Mu) +
      7200*Quad(g1d)*Sqr(Mu) + 1600*Quad(g1dp)*Sqr(Mu) - 7200*Quad(g2)*Sqr(Mu)
      - 720*Sqr(g1)*Sqr(g1d)*Sqr(Mu) - 240*Sqr(g1)*Sqr(g1dp)*Sqr(Mu) + 2400*Sqr
      (g1d)*Sqr(g1dp)*Sqr(Mu) - 8400*Sqr(g1d)*Sqr(g2)*Sqr(Mu) - 1200*Sqr(g1dp)*
      Sqr(g2)*Sqr(Mu) + 4200*Sqr(g1d)*Sqr(g2u)*Sqr(Mu) + 600*Sqr(g1dp)*Sqr(g2u)
      *Sqr(Mu) + 600*Sqr(g1d)*Sqr(g2up)*Sqr(Mu) + 1800*Sqr(g1dp)*Sqr(g2up)*Sqr(
      Mu)));
   const double beta_M112_2 = Re(0.5*twoLoop*(9*Conj(Lambda5)*Conj(Lambda7)*
      Conj(M122) + 12*g1dp*g2up*Lambda6*Conj(MassB)*Mu + 36*g1d*g2u*Lambda6*
      Conj(MassWB)*Mu + 39*AbsSqr(MassWB)*Quad(g1d) + 11*AbsSqr(MassB)*Quad(
      g1dp) - 72*AbsSqr(MassWB)*Quad(g2) + 9*AbsSqr(MassB)*Sqr(g1d)*Sqr(g1dp) +
      9*AbsSqr(MassWB)*Sqr(g1d)*Sqr(g1dp) + 6*MassWB*Conj(MassB)*Sqr(g1d)*Sqr(
      g1dp) + 6*MassB*Conj(MassWB)*Sqr(g1d)*Sqr(g1dp) - 72*AbsSqr(MassWB)*Sqr(
      g1d)*Sqr(g2) + 18*AbsSqr(MassWB)*Sqr(g1d)*Sqr(g2u) + 6*AbsSqr(MassB)*Sqr(
      g1dp)*Sqr(g2up)));

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

/**
 * Calculates the 5-loop beta function of M112.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M112_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M112;

   beta_M112 = 0;


   return beta_M112;
}

} // namespace flexiblesusy
