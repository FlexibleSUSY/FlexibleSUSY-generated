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

// File generated at Wed 16 Oct 2019 22:00:05

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

   beta_M122 = Re(0.1*oneOver16PiSqr*(-60*Lambda6*M112 + 20*Lambda3*M122 + 40*
      Lambda4*M122 - 60*Lambda7*M222 + 30*M122*traceYdAdjYd + 10*M122*
      traceYeAdjYe + 30*M122*traceYuAdjYu + 60*Conj(Lambda5)*Conj(M122) + 20*
      g1dp*g2up*MassB*Mu + 60*g1d*g2u*MassWB*Mu - 9*M122*Sqr(g1) + 15*M122*Sqr(
      g1d) + 5*M122*Sqr(g1dp) - 45*M122*Sqr(g2) + 15*M122*Sqr(g2u) + 5*M122*Sqr
      (g2up)));


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

   const double beta_M122_1 = Re(0.0025*twoLoop*(13200*Lambda1*Lambda6*M112 +
      4200*Lambda3*Lambda6*M112 + 4200*Lambda4*Lambda6*M112 - 1200*Lambda2*
      Lambda7*M112 + 1800*Lambda3*Lambda7*M112 + 1800*Lambda4*Lambda7*M112 +
      2400*g1d*g1dp*g2u*g2up*M122 - 4800*Lambda1*Lambda3*M122 - 4800*Lambda2*
      Lambda3*M122 - 4800*Lambda1*Lambda4*M122 - 4800*Lambda2*Lambda4*M122 -
      2400*Lambda3*Lambda4*M122 - 1200*Lambda1*Lambda6*M222 + 1800*Lambda3*
      Lambda6*M222 + 1800*Lambda4*Lambda6*M222 + 13200*Lambda2*Lambda7*M222 +
      4200*Lambda3*Lambda7*M222 + 4200*Lambda4*Lambda7*M222 + 14400*Lambda6*
      M112*traceYdAdjYd - 2400*Lambda3*M122*traceYdAdjYd - 4800*Lambda4*M122*
      traceYdAdjYd - 2700*M122*traceYdAdjYdYdAdjYd - 6600*M122*
      traceYdAdjYuYuAdjYd + 4800*Lambda6*M112*traceYeAdjYe - 800*Lambda3*M122*
      traceYeAdjYe - 1600*Lambda4*M122*traceYeAdjYe - 900*M122*
      traceYeAdjYeYeAdjYe - 2400*Lambda3*M122*traceYuAdjYu - 4800*Lambda4*M122*
      traceYuAdjYu + 14400*Lambda7*M222*traceYuAdjYu - 2700*M122*
      traceYuAdjYuYuAdjYu - 800*g1dp*g2up*Lambda3*MassB*Mu - 1600*g1dp*g2up*
      Lambda4*MassB*Mu - 2400*g1d*g2u*Lambda3*MassWB*Mu - 4800*g1d*g2u*Lambda4*
      MassWB*Mu - 3000*g2u*MassWB*Cube(g1d)*Mu - 1400*g2up*MassB*Cube(g1dp)*Mu
      - 3000*g1d*MassWB*Cube(g2u)*Mu - 1400*g1dp*MassB*Cube(g2up)*Mu + 1497*
      M122*Quad(g1) - 1125*M122*Quad(g1d) - 225*M122*Quad(g1dp) - 3075*M122*
      Quad(g2) - 1125*M122*Quad(g2u) - 225*M122*Quad(g2up) - 2880*Lambda6*M112*
      Sqr(g1) + 960*Lambda3*M122*Sqr(g1) + 1920*Lambda4*M122*Sqr(g1) - 2880*
      Lambda7*M222*Sqr(g1) + 250*M122*traceYdAdjYd*Sqr(g1) + 750*M122*
      traceYeAdjYe*Sqr(g1) + 850*M122*traceYuAdjYu*Sqr(g1) + 120*g1dp*g2up*
      MassB*Mu*Sqr(g1) + 360*g1d*g2u*MassWB*Mu*Sqr(g1) + 7200*Lambda6*M112*Sqr(
      g1d) - 1200*Lambda3*M122*Sqr(g1d) - 2400*Lambda4*M122*Sqr(g1d) - 600*g1dp
      *g2up*MassB*Mu*Sqr(g1d) - 1200*g1dp*g2up*MassWB*Mu*Sqr(g1d) + 225*M122*
      Sqr(g1)*Sqr(g1d) + 2400*Lambda6*M112*Sqr(g1dp) - 400*Lambda3*M122*Sqr(
      g1dp) - 800*Lambda4*M122*Sqr(g1dp) - 1200*g1d*g2u*MassB*Mu*Sqr(g1dp) -
      600*g1d*g2u*MassWB*Mu*Sqr(g1dp) + 75*M122*Sqr(g1)*Sqr(g1dp) - 450*M122*
      Sqr(g1d)*Sqr(g1dp) - 14400*Lambda6*M112*Sqr(g2) + 4800*Lambda3*M122*Sqr(
      g2) + 9600*Lambda4*M122*Sqr(g2) - 14400*Lambda7*M222*Sqr(g2) + 2250*M122*
      traceYdAdjYd*Sqr(g2) + 750*M122*traceYeAdjYe*Sqr(g2) + 2250*M122*
      traceYuAdjYu*Sqr(g2) + 600*g1dp*g2up*MassB*Mu*Sqr(g2) + 11400*g1d*g2u*
      MassWB*Mu*Sqr(g2) + 450*M122*Sqr(g1)*Sqr(g2) + 4125*M122*Sqr(g1d)*Sqr(g2)
      + 375*M122*Sqr(g1dp)*Sqr(g2) - 1200*Lambda3*M122*Sqr(g2u) - 2400*Lambda4*
      M122*Sqr(g2u) + 7200*Lambda7*M222*Sqr(g2u) - 600*g1dp*g2up*MassB*Mu*Sqr(
      g2u) - 1200*g1dp*g2up*MassWB*Mu*Sqr(g2u) + 225*M122*Sqr(g1)*Sqr(g2u) -
      2100*M122*Sqr(g1d)*Sqr(g2u) + 4125*M122*Sqr(g2)*Sqr(g2u) - 400*Lambda3*
      M122*Sqr(g2up) - 800*Lambda4*M122*Sqr(g2up) + 2400*Lambda7*M222*Sqr(g2up)
      - 1200*g1d*g2u*MassB*Mu*Sqr(g2up) - 600*g1d*g2u*MassWB*Mu*Sqr(g2up) + 75*
      M122*Sqr(g1)*Sqr(g2up) + 100*M122*Sqr(g1dp)*Sqr(g2up) + 375*M122*Sqr(g2)*
      Sqr(g2up) - 450*M122*Sqr(g2u)*Sqr(g2up) + 8000*M122*traceYdAdjYd*Sqr(g3)
      + 8000*M122*traceYuAdjYu*Sqr(g3) + 2400*M122*Sqr(Lambda1) + 2400*M122*Sqr
      (Lambda2)));
   const double beta_M122_2 = Re(0.3*twoLoop*(10*M122*AbsSqr(Lambda5) - 40*
      Lambda7*M122*Conj(Lambda6) + 35*M112*Conj(Lambda5)*Conj(Lambda6) + 15*
      M222*Conj(Lambda5)*Conj(Lambda6) - 40*Lambda6*M122*Conj(Lambda7) + 15*
      M112*Conj(Lambda5)*Conj(Lambda7) + 35*M222*Conj(Lambda5)*Conj(Lambda7) -
      40*Lambda6*Lambda7*Conj(M122) - 40*Lambda1*Conj(Lambda5)*Conj(M122) - 40*
      Lambda2*Conj(Lambda5)*Conj(M122) - 40*Lambda3*Conj(Lambda5)*Conj(M122) -
      40*Lambda4*Conj(Lambda5)*Conj(M122) - 60*traceYdAdjYd*Conj(Lambda5)*Conj(
      M122) - 20*traceYeAdjYe*Conj(Lambda5)*Conj(M122) - 60*traceYuAdjYu*Conj(
      Lambda5)*Conj(M122) - 20*g1dp*g2up*Conj(Lambda5)*Conj(MassB)*Mu - 60*g1d*
      g2u*Conj(Lambda5)*Conj(MassWB)*Mu + 24*Conj(Lambda5)*Conj(M122)*Sqr(g1) -
      30*Conj(Lambda5)*Conj(M122)*Sqr(g1d) - 10*Conj(Lambda5)*Conj(M122)*Sqr(
      g1dp) + 120*Conj(Lambda5)*Conj(M122)*Sqr(g2) - 30*Conj(Lambda5)*Conj(M122
      )*Sqr(g2u) - 10*Conj(Lambda5)*Conj(M122)*Sqr(g2up) - 40*Conj(M122)*Sqr(
      Lambda6) - 40*Conj(M122)*Sqr(Lambda7)));

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

/**
 * Calculates the 5-loop beta function of M122.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M122_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M122;

   beta_M122 = 0;


   return beta_M122;
}

} // namespace flexiblesusy
