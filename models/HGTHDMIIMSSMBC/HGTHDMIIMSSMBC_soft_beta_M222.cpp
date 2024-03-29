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
 * Calculates the 1-loop beta function of M222.
 *
 * @return 1-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M222_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M222;

   beta_M222 = Re(0.1*(40*Lambda3*M112 + 20*Lambda4*M112 + 120*Lambda2*M222 +
      60*M222*traceYuAdjYu - 60*M122*Conj(Lambda7) - 60*Lambda7*Conj(M122) - 9*
      M222*Sqr(g1) - 45*M222*Sqr(g2) + 30*M222*Sqr(g2u) - 60*AbsSqr(MassWB)*Sqr
      (g2u) + 10*M222*Sqr(g2up) - 20*AbsSqr(MassB)*Sqr(g2up) - 60*Sqr(g2u)*Sqr(
      Mu) - 20*Sqr(g2up)*Sqr(Mu)));


   return oneLoop * beta_M222;
}

/**
 * Calculates the 2-loop beta function of M222.
 *
 * @return 2-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M222_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M222;

   const double beta_M222_1 = Re(0.0025*(-3200*Lambda3*Lambda4*M112 + 1800*
      Lambda5*Lambda6*M122 + 4200*Lambda5*Lambda7*M122 - 800*Lambda3*Lambda4*
      M222 - 9600*Lambda3*M112*traceYdAdjYd - 4800*Lambda4*M112*traceYdAdjYd -
      1800*M222*traceYdAdjYuYuAdjYd - 3200*Lambda3*M112*traceYeAdjYe - 1600*
      Lambda4*M112*traceYeAdjYe - 28800*Lambda2*M222*traceYuAdjYu - 5400*M222*
      traceYuAdjYuYuAdjYu - 4800*M112*AbsSqr(Lambda5) - 1200*M222*AbsSqr(
      Lambda5) - 7200*M112*AbsSqr(Lambda6) + 1200*M222*AbsSqr(Lambda6) - 7200*
      M112*AbsSqr(Lambda7) - 10800*M222*AbsSqr(Lambda7) - 1200*Lambda1*M122*
      Conj(Lambda6) + 1800*Lambda3*M122*Conj(Lambda6) + 1800*Lambda4*M122*Conj(
      Lambda6) + 13200*Lambda2*M122*Conj(Lambda7) + 4200*Lambda3*M122*Conj(
      Lambda7) + 4200*Lambda4*M122*Conj(Lambda7) + 7200*M122*traceYdAdjYd*Conj(
      Lambda7) + 2400*M122*traceYeAdjYe*Conj(Lambda7) + 7200*M122*traceYuAdjYu*
      Conj(Lambda7) - 1200*Lambda1*Lambda6*Conj(M122) + 1800*Lambda3*Lambda6*
      Conj(M122) + 1800*Lambda4*Lambda6*Conj(M122) + 13200*Lambda2*Lambda7*Conj
      (M122) + 4200*Lambda3*Lambda7*Conj(M122) + 4200*Lambda4*Lambda7*Conj(M122
      ) + 7200*Lambda7*traceYdAdjYd*Conj(M122) + 2400*Lambda7*traceYeAdjYe*Conj
      (M122) + 7200*Lambda7*traceYuAdjYu*Conj(M122) + 1800*Conj(Lambda5)*Conj(
      Lambda6)*Conj(M122) + 4200*Conj(Lambda5)*Conj(Lambda7)*Conj(M122) + 2400*
      g1dp*g2up*MassB*Conj(Lambda7)*Mu + 7200*g1d*g2u*MassWB*Conj(Lambda7)*Mu +
      360*M112*Quad(g1) + 1857*M222*Quad(g1) + 3000*M112*Quad(g2) - 75*M222*
      Quad(g2) - 2250*M222*Quad(g2u) - 450*M222*Quad(g2up) + 1920*Lambda3*M112*
      Sqr(g1) + 960*Lambda4*M112*Sqr(g1) + 5760*Lambda2*M222*Sqr(g1) + 1700*
      M222*traceYuAdjYu*Sqr(g1) - 2880*M122*Conj(Lambda7)*Sqr(g1) - 2880*
      Lambda7*Conj(M122)*Sqr(g1) - 4800*Lambda3*M112*Sqr(g1d) - 2400*Lambda4*
      M112*Sqr(g1d) + 3600*M122*Conj(Lambda7)*Sqr(g1d) + 3600*Lambda7*Conj(M122
      )*Sqr(g1d) - 1600*Lambda3*M112*Sqr(g1dp) - 800*Lambda4*M112*Sqr(g1dp) +
      1200*M122*Conj(Lambda7)*Sqr(g1dp) + 1200*Lambda7*Conj(M122)*Sqr(g1dp) +
      9600*Lambda3*M112*Sqr(g2) + 4800*Lambda4*M112*Sqr(g2) + 28800*Lambda2*
      M222*Sqr(g2) + 4500*M222*traceYuAdjYu*Sqr(g2) - 14400*M122*Conj(Lambda7)*
      Sqr(g2) - 14400*Lambda7*Conj(M122)*Sqr(g2) + 450*M222*Sqr(g1)*Sqr(g2) -
      14400*Lambda2*M222*Sqr(g2u) + 3600*M122*Conj(Lambda7)*Sqr(g2u) + 3600*
      Lambda7*Conj(M122)*Sqr(g2u) + 450*M222*Sqr(g1)*Sqr(g2u) - 900*M222*Sqr(
      g1d)*Sqr(g2u) + 8250*M222*Sqr(g2)*Sqr(g2u) - 4800*Lambda2*M222*Sqr(g2up)
      + 1200*M122*Conj(Lambda7)*Sqr(g2up) + 1200*Lambda7*Conj(M122)*Sqr(g2up) +
      150*M222*Sqr(g1)*Sqr(g2up) - 300*M222*Sqr(g1dp)*Sqr(g2up) + 1200*AbsSqr(
      MassB)*Sqr(g1dp)*Sqr(g2up) + 750*M222*Sqr(g2)*Sqr(g2up) - 900*M222*Sqr(
      g2u)*Sqr(g2up) + 16000*M222*traceYuAdjYu*Sqr(g3) - 24000*M222*Sqr(Lambda2
      ) - 3200*M112*Sqr(Lambda3) - 800*M222*Sqr(Lambda3) - 3200*M112*Sqr(
      Lambda4) - 800*M222*Sqr(Lambda4) + 2400*g1d*g1dp*g2u*g2up*Sqr(Mu) - 864*
      Quad(g1)*Sqr(Mu) - 7200*Quad(g2)*Sqr(Mu) + 7200*Quad(g2u)*Sqr(Mu) + 1600*
      Quad(g2up)*Sqr(Mu) - 720*Sqr(g1)*Sqr(g2u)*Sqr(Mu) + 4200*Sqr(g1d)*Sqr(g2u
      )*Sqr(Mu) + 600*Sqr(g1dp)*Sqr(g2u)*Sqr(Mu) - 8400*Sqr(g2)*Sqr(g2u)*Sqr(Mu
      ) - 240*Sqr(g1)*Sqr(g2up)*Sqr(Mu) + 600*Sqr(g1d)*Sqr(g2up)*Sqr(Mu) + 1800
      *Sqr(g1dp)*Sqr(g2up)*Sqr(Mu) - 1200*Sqr(g2)*Sqr(g2up)*Sqr(Mu) + 2400*Sqr(
      g2u)*Sqr(g2up)*Sqr(Mu)));
   const double beta_M222_2 = Re(0.5*(12*g1dp*g2up*Lambda7*Conj(MassB)*Mu + 36*
      g1d*g2u*Lambda7*Conj(MassWB)*Mu - 72*AbsSqr(MassWB)*Quad(g2) + 39*AbsSqr(
      MassWB)*Quad(g2u) + 11*AbsSqr(MassB)*Quad(g2up) + 18*AbsSqr(MassWB)*Sqr(
      g1d)*Sqr(g2u) - 72*AbsSqr(MassWB)*Sqr(g2)*Sqr(g2u) + 9*AbsSqr(MassB)*Sqr(
      g2u)*Sqr(g2up) + 9*AbsSqr(MassWB)*Sqr(g2u)*Sqr(g2up) + 6*MassWB*Conj(
      MassB)*Sqr(g2u)*Sqr(g2up) + 6*MassB*Conj(MassWB)*Sqr(g2u)*Sqr(g2up)));

   beta_M222 = beta_M222_1 + beta_M222_2;


   return twoLoop * beta_M222;
}

/**
 * Calculates the 3-loop beta function of M222.
 *
 * @return 3-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M222_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return threeLoop * beta_M222;
}

/**
 * Calculates the 4-loop beta function of M222.
 *
 * @return 4-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M222_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return fourLoop * beta_M222;
}

/**
 * Calculates the 5-loop beta function of M222.
 *
 * @return 5-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M222_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M222;

   beta_M222 = 0;


   return fiveLoop * beta_M222;
}

} // namespace flexiblesusy
