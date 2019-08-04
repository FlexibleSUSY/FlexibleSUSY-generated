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

// File generated at Sun 4 Aug 2019 19:24:39

#include "SplitMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mu2.
 *
 * @return 1-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_mu2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_mu2;

   beta_mu2 = Re(0.1*oneOver16PiSqr*(60*mu2*traceYdAdjYd + 20*mu2*traceYeAdjYe
      + 60*mu2*traceYuAdjYu + 60*mu2*Lambdax + 40*gYd*gYu*MassB*Mu + 120*g2d*
      g2u*MassWB*Mu - 9*mu2*Sqr(g1) - 45*mu2*Sqr(g2) + 30*mu2*Sqr(g2d) + 30*mu2
      *Sqr(g2u) + 10*mu2*Sqr(gYd) + 10*mu2*Sqr(gYu) + 20*Sqr(gYd)*Sqr(MassB) +
      20*Sqr(gYu)*Sqr(MassB) + 60*Sqr(g2d)*Sqr(MassWB) + 60*Sqr(g2u)*Sqr(MassWB
      ) + 60*Sqr(g2d)*Sqr(Mu) + 60*Sqr(g2u)*Sqr(Mu) + 20*Sqr(gYd)*Sqr(Mu) + 20*
      Sqr(gYu)*Sqr(Mu)));


   return beta_mu2;
}

/**
 * Calculates the 2-loop beta function of mu2.
 *
 * @return 2-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_mu2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_mu2;

   beta_mu2 = Re(0.0025*twoLoop*(-4800*g2d*g2u*gYd*gYu*MassB*MassWB + 2400*g2d*
      g2u*gYd*gYu*mu2 - 5400*mu2*traceYdAdjYdYdAdjYd - 8400*mu2*
      traceYdAdjYuYuAdjYd - 1800*mu2*traceYeAdjYeYeAdjYe - 5400*mu2*
      traceYuAdjYuYuAdjYu - 14400*mu2*traceYdAdjYd*Lambdax - 4800*mu2*
      traceYeAdjYe*Lambdax - 14400*mu2*traceYuAdjYu*Lambdax - 20400*g2u*MassWB*
      Cube(g2d)*Mu - 20400*g2d*MassWB*Cube(g2u)*Mu - 7600*gYu*MassB*Cube(gYd)*
      Mu - 7600*gYd*MassB*Cube(gYu)*Mu - 4800*gYd*gYu*MassB*Lambdax*Mu - 14400*
      g2d*g2u*MassWB*Lambdax*Mu + 1791*mu2*Quad(g1) - 625*mu2*Quad(g2) - 2250*
      mu2*Quad(g2d) - 2250*mu2*Quad(g2u) - 450*mu2*Quad(gYd) - 450*mu2*Quad(gYu
      ) + 500*mu2*traceYdAdjYd*Sqr(g1) + 1500*mu2*traceYeAdjYe*Sqr(g1) + 1700*
      mu2*traceYuAdjYu*Sqr(g1) + 2880*mu2*Lambdax*Sqr(g1) + 240*gYd*gYu*MassB*
      Mu*Sqr(g1) + 720*g2d*g2u*MassWB*Mu*Sqr(g1) + 4500*mu2*traceYdAdjYd*Sqr(g2
      ) + 1500*mu2*traceYeAdjYe*Sqr(g2) + 4500*mu2*traceYuAdjYu*Sqr(g2) + 14400
      *mu2*Lambdax*Sqr(g2) + 1200*gYd*gYu*MassB*Mu*Sqr(g2) + 22800*g2d*g2u*
      MassWB*Mu*Sqr(g2) + 450*mu2*Sqr(g1)*Sqr(g2) - 7200*mu2*Lambdax*Sqr(g2d) -
      3600*gYd*gYu*MassB*Mu*Sqr(g2d) - 4800*gYd*gYu*MassWB*Mu*Sqr(g2d) + 450*
      mu2*Sqr(g1)*Sqr(g2d) + 8250*mu2*Sqr(g2)*Sqr(g2d) - 7200*mu2*Lambdax*Sqr(
      g2u) - 3600*gYd*gYu*MassB*Mu*Sqr(g2u) - 4800*gYd*gYu*MassWB*Mu*Sqr(g2u) +
      450*mu2*Sqr(g1)*Sqr(g2u) + 8250*mu2*Sqr(g2)*Sqr(g2u) - 3000*mu2*Sqr(g2d)*
      Sqr(g2u) + 16000*mu2*traceYdAdjYd*Sqr(g3) + 16000*mu2*traceYuAdjYu*Sqr(g3
      ) - 2400*mu2*Lambdax*Sqr(gYd) - 4800*g2d*g2u*MassB*Mu*Sqr(gYd) - 3600*g2d
      *g2u*MassWB*Mu*Sqr(gYd) + 150*mu2*Sqr(g1)*Sqr(gYd) + 750*mu2*Sqr(g2)*Sqr(
      gYd) - 2400*MassB*MassWB*Sqr(g2d)*Sqr(gYd) - 900*mu2*Sqr(g2d)*Sqr(gYd) -
      2400*mu2*Lambdax*Sqr(gYu) - 4800*g2d*g2u*MassB*Mu*Sqr(gYu) - 3600*g2d*g2u
      *MassWB*Mu*Sqr(gYu) + 150*mu2*Sqr(g1)*Sqr(gYu) + 750*mu2*Sqr(g2)*Sqr(gYu)
      - 2400*MassB*MassWB*Sqr(g2u)*Sqr(gYu) - 900*mu2*Sqr(g2u)*Sqr(gYu) - 200*
      mu2*Sqr(gYd)*Sqr(gYu) - 4800*g2d*g2u*gYd*gYu*Sqr(MassB) - 2200*Quad(gYd)*
      Sqr(MassB) - 2200*Quad(gYu)*Sqr(MassB) - 1800*Sqr(g2d)*Sqr(gYd)*Sqr(MassB
      ) - 1800*Sqr(g2u)*Sqr(gYu)*Sqr(MassB) - 4800*Sqr(gYd)*Sqr(gYu)*Sqr(MassB)
      - 4800*g2d*g2u*gYd*gYu*Sqr(MassWB) + 14400*Quad(g2)*Sqr(MassWB) - 7800*
      Quad(g2d)*Sqr(MassWB) - 7800*Quad(g2u)*Sqr(MassWB) + 14400*Sqr(g2)*Sqr(
      g2d)*Sqr(MassWB) + 14400*Sqr(g2)*Sqr(g2u)*Sqr(MassWB) - 9600*Sqr(g2d)*Sqr
      (g2u)*Sqr(MassWB) - 1800*Sqr(g2d)*Sqr(gYd)*Sqr(MassWB) - 1800*Sqr(g2u)*
      Sqr(gYu)*Sqr(MassWB) - 6000*mu2*Sqr(Lambdax) - 14400*g2d*g2u*gYd*gYu*Sqr(
      Mu) + 864*Quad(g1)*Sqr(Mu) + 7200*Quad(g2)*Sqr(Mu) - 7200*Quad(g2d)*Sqr(
      Mu) - 7200*Quad(g2u)*Sqr(Mu) - 1600*Quad(gYd)*Sqr(Mu) - 1600*Quad(gYu)*
      Sqr(Mu) + 720*Sqr(g1)*Sqr(g2d)*Sqr(Mu) + 8400*Sqr(g2)*Sqr(g2d)*Sqr(Mu) +
      720*Sqr(g1)*Sqr(g2u)*Sqr(Mu) + 8400*Sqr(g2)*Sqr(g2u)*Sqr(Mu) - 13200*Sqr(
      g2d)*Sqr(g2u)*Sqr(Mu) + 240*Sqr(g1)*Sqr(gYd)*Sqr(Mu) + 1200*Sqr(g2)*Sqr(
      gYd)*Sqr(Mu) - 2400*Sqr(g2d)*Sqr(gYd)*Sqr(Mu) - 1200*Sqr(g2u)*Sqr(gYd)*
      Sqr(Mu) + 240*Sqr(g1)*Sqr(gYu)*Sqr(Mu) + 1200*Sqr(g2)*Sqr(gYu)*Sqr(Mu) -
      1200*Sqr(g2d)*Sqr(gYu)*Sqr(Mu) - 2400*Sqr(g2u)*Sqr(gYu)*Sqr(Mu) - 8400*
      Sqr(gYd)*Sqr(gYu)*Sqr(Mu)));


   return beta_mu2;
}

/**
 * Calculates the 3-loop beta function of mu2.
 *
 * @return 3-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_mu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mu2;

   beta_mu2 = 0;


   return beta_mu2;
}

/**
 * Calculates the 4-loop beta function of mu2.
 *
 * @return 4-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_mu2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mu2;

   beta_mu2 = 0;


   return beta_mu2;
}

/**
 * Calculates the 5-loop beta function of mu2.
 *
 * @return 5-loop beta function
 */
double SplitMSSM_soft_parameters::calc_beta_mu2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mu2;

   beta_mu2 = 0;


   return beta_mu2;
}

} // namespace flexiblesusy
