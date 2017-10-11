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

// File generated at Tue 10 Oct 2017 21:10:03

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

   beta_mu2 = Re(oneOver16PiSqr*(6*mu2*traceYdAdjYd + 2*mu2*traceYeAdjYe
      + 6*mu2*traceYuAdjYu + 6*mu2*Lambdax + 2*gYd*gYu*MassB*Mu + 6*g2d*g2u*
      MassWB*Mu - 0.9*mu2*Sqr(g1) - 4.5*mu2*Sqr(g2) + 3*mu2*Sqr(g2d) + 3*mu2*
      Sqr(g2u) + 6*Conj(MassWB)*(g2d*g2u*Mu + MassWB*Sqr(g2d) + MassWB*Sqr(g2u)
      ) + mu2*Sqr(gYd) + mu2*Sqr(gYu) + 2*Conj(MassB)*(gYd*gYu*Mu + MassB*Sqr(
      gYd) + MassB*Sqr(gYu)) + 6*Sqr(g2d)*Sqr(Mu) + 6*Sqr(g2u)*Sqr(Mu) + 2*Sqr(
      gYd)*Sqr(Mu) + 2*Sqr(gYu)*Sqr(Mu)));


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

   const double beta_mu2_1 = Re(twoLoop*(6*g2d*g2u*gYd*gYu*mu2 - 13.5*mu2
      *traceYdAdjYdYdAdjYd - 21*mu2*traceYdAdjYuYuAdjYd - 4.5*mu2*
      traceYeAdjYeYeAdjYe - 13.5*mu2*traceYuAdjYuYuAdjYu - 36*mu2*traceYdAdjYd*
      Lambdax - 12*mu2*traceYeAdjYe*Lambdax - 36*mu2*traceYuAdjYu*Lambdax -
      25.5*g2u*MassWB*Cube(g2d)*Mu - 25.5*g2d*MassWB*Cube(g2u)*Mu - 9.5*gYu*
      MassB*Cube(gYd)*Mu - 9.5*gYd*MassB*Cube(gYu)*Mu - 6*gYd*gYu*MassB*Lambdax
      *Mu - 18*g2d*g2u*MassWB*Lambdax*Mu + 4.4775*mu2*Quad(g1) - 1.5625*mu2*
      Quad(g2) - 5.625*mu2*Quad(g2d) - 5.625*mu2*Quad(g2u) - 1.125*mu2*Quad(gYd
      ) - 1.125*mu2*Quad(gYu) + 1.25*mu2*traceYdAdjYd*Sqr(g1) + 3.75*mu2*
      traceYeAdjYe*Sqr(g1) + 4.25*mu2*traceYuAdjYu*Sqr(g1) + 7.2*mu2*Lambdax*
      Sqr(g1) + 0.3*gYd*gYu*MassB*Mu*Sqr(g1) + 0.9*g2d*g2u*MassWB*Mu*Sqr(g1) +
      11.25*mu2*traceYdAdjYd*Sqr(g2) + 3.75*mu2*traceYeAdjYe*Sqr(g2) + 11.25*
      mu2*traceYuAdjYu*Sqr(g2) + 36*mu2*Lambdax*Sqr(g2) + 1.5*gYd*gYu*MassB*Mu*
      Sqr(g2) + 28.5*g2d*g2u*MassWB*Mu*Sqr(g2) + 1.125*mu2*Sqr(g1)*Sqr(g2) - 18
      *mu2*Lambdax*Sqr(g2d) - 4.5*gYd*gYu*MassB*Mu*Sqr(g2d) - 6*gYd*gYu*MassWB*
      Mu*Sqr(g2d) + 1.125*mu2*Sqr(g1)*Sqr(g2d) + 20.625*mu2*Sqr(g2)*Sqr(g2d) -
      18*mu2*Lambdax*Sqr(g2u) - 4.5*gYd*gYu*MassB*Mu*Sqr(g2u) - 6*gYd*gYu*
      MassWB*Mu*Sqr(g2u) + 1.125*mu2*Sqr(g1)*Sqr(g2u) + 20.625*mu2*Sqr(g2)*Sqr(
      g2u) - 7.5*mu2*Sqr(g2d)*Sqr(g2u) + 40*mu2*traceYdAdjYd*Sqr(g3) + 40*mu2*
      traceYuAdjYu*Sqr(g3) - 6*mu2*Lambdax*Sqr(gYd) - 6*g2d*g2u*MassB*Mu*Sqr(
      gYd) - 4.5*g2d*g2u*MassWB*Mu*Sqr(gYd) + 0.375*mu2*Sqr(g1)*Sqr(gYd) +
      1.875*mu2*Sqr(g2)*Sqr(gYd) - 2.25*mu2*Sqr(g2d)*Sqr(gYd) - 3*MassB*Conj(
      MassWB)*Sqr(g2d)*Sqr(gYd) - 6*mu2*Lambdax*Sqr(gYu) - 6*g2d*g2u*MassB*Mu*
      Sqr(gYu) - 4.5*g2d*g2u*MassWB*Mu*Sqr(gYu) + 0.375*mu2*Sqr(g1)*Sqr(gYu) +
      1.875*mu2*Sqr(g2)*Sqr(gYu) - 2.25*mu2*Sqr(g2u)*Sqr(gYu) - 0.5*mu2*Sqr(gYd
      )*Sqr(gYu) - 0.1*Conj(MassB)*(95*gYu*Cube(gYd)*Mu + 55*MassB*Quad(gYd) +
      15*gYd*(3*gYd*MassB + 2*gYd*MassWB + 3*gYu*Mu)*Sqr(g2d) + 120*MassB*Sqr(
      gYd)*Sqr(gYu) + gYd*gYu*Mu*(60*Lambdax - 3*Sqr(g1) - 15*Sqr(g2) + 45*Sqr(
      g2u) + 95*Sqr(gYu)) + 5*Sqr(gYu)*((9*MassB + 6*MassWB)*Sqr(g2u) + 11*
      MassB*Sqr(gYu)) + 60*g2d*g2u*(gYd*gYu*(2*MassB + MassWB) + Mu*Sqr(gYd) +
      Mu*Sqr(gYu))) - 15*mu2*Sqr(Lambdax) - 36*g2d*g2u*gYd*gYu*Sqr(Mu) + 2.16*
      Quad(g1)*Sqr(Mu) + 18*Quad(g2)*Sqr(Mu) - 18*Quad(g2d)*Sqr(Mu) - 18*Quad(
      g2u)*Sqr(Mu) - 4*Quad(gYd)*Sqr(Mu) - 4*Quad(gYu)*Sqr(Mu) + 1.8*Sqr(g1)*
      Sqr(g2d)*Sqr(Mu) + 21*Sqr(g2)*Sqr(g2d)*Sqr(Mu) + 1.8*Sqr(g1)*Sqr(g2u)*Sqr
      (Mu) + 21*Sqr(g2)*Sqr(g2u)*Sqr(Mu) - 33*Sqr(g2d)*Sqr(g2u)*Sqr(Mu) + 0.6*
      Sqr(g1)*Sqr(gYd)*Sqr(Mu) + 3*Sqr(g2)*Sqr(gYd)*Sqr(Mu) - 6*Sqr(g2d)*Sqr(
      gYd)*Sqr(Mu) - 3*Sqr(g2u)*Sqr(gYd)*Sqr(Mu) + 0.6*Sqr(g1)*Sqr(gYu)*Sqr(Mu)
      + 3*Sqr(g2)*Sqr(gYu)*Sqr(Mu) - 3*Sqr(g2d)*Sqr(gYu)*Sqr(Mu) - 6*Sqr(g2u)*
      Sqr(gYu)*Sqr(Mu) - 21*Sqr(gYd)*Sqr(gYu)*Sqr(Mu)));
   const double beta_mu2_2 = Re(-0.3*twoLoop*Conj(MassWB)*(85*g2u*Cube(
      g2d)*Mu + 65*MassWB*Quad(g2d) + 5*Sqr(g2d)*(4*gYd*gYu*Mu - 24*MassWB*Sqr(
      g2) + 16*MassWB*Sqr(g2u) + 3*MassWB*Sqr(gYd)) + g2d*g2u*(20*gYd*gYu*(
      MassB + 2*MassWB) + 15*Mu*Sqr(gYd) + Mu*(60*Lambdax - 3*Sqr(g1) - 95*Sqr(
      g2) + 85*Sqr(g2u) + 15*Sqr(gYu))) + 5*(-24*MassWB*Quad(g2) + 13*MassWB*
      Quad(g2u) + Sqr(g2u)*(4*gYd*gYu*Mu - 24*MassWB*Sqr(g2) + (2*MassB + 3*
      MassWB)*Sqr(gYu)))));

   beta_mu2 = beta_mu2_1 + beta_mu2_2;


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

} // namespace flexiblesusy