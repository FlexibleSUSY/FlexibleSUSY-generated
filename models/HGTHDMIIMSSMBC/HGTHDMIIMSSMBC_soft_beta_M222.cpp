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

// File generated at Tue 10 Oct 2017 20:54:05

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

   beta_M222 = Re(oneOver16PiSqr*(4*Lambda3*M112 + 2*Lambda4*M112 - 12*
      Lambda7*M122 + 12*Lambda2*M222 + 6*M222*traceYuAdjYu - 0.9*M222*Sqr(g1) -
      4.5*M222*Sqr(g2) + 3*M222*Sqr(g2u) - 6*AbsSqr(MassWB)*Sqr(g2u) + M222*
      Sqr(g2up) - 2*AbsSqr(MassB)*Sqr(g2up) - 6*Sqr(g2u)*Sqr(Mu) - 2*Sqr(g2up)*
      Sqr(Mu)));


   return beta_M222;
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

   beta_M222 = Re(twoLoop*(-8*Lambda3*Lambda4*M112 + 12*Lambda3*Lambda6*
      M122 + 12*Lambda4*Lambda6*M122 + 12*Lambda5*Lambda6*M122 + 72*Lambda2*
      Lambda7*M122 + 24*Lambda3*Lambda7*M122 + 24*Lambda4*Lambda7*M122 + 24*
      Lambda5*Lambda7*M122 - 2*Lambda3*Lambda4*M222 - 24*Lambda3*M112*
      traceYdAdjYd - 12*Lambda4*M112*traceYdAdjYd + 36*Lambda7*M122*
      traceYdAdjYd - 4.5*M222*traceYdAdjYuYuAdjYd - 8*Lambda3*M112*traceYeAdjYe
      - 4*Lambda4*M112*traceYeAdjYe + 12*Lambda7*M122*traceYeAdjYe + 36*
      Lambda7*M122*traceYuAdjYu - 72*Lambda2*M222*traceYuAdjYu - 13.5*M222*
      traceYuAdjYuYuAdjYu + 6*g1dp*g2up*Lambda7*MassB*Mu + 18*g1d*g2u*Lambda7*
      MassWB*Mu + 0.9*M112*Quad(g1) + 4.6425*M222*Quad(g1) + 7.5*M112*Quad(g2)
      - 0.1875*M222*Quad(g2) - 5.625*M222*Quad(g2u) - 1.125*M222*Quad(g2up) +
      4.8*Lambda3*M112*Sqr(g1) + 2.4*Lambda4*M112*Sqr(g1) - 14.4*Lambda7*M122*
      Sqr(g1) + 14.4*Lambda2*M222*Sqr(g1) + 4.25*M222*traceYuAdjYu*Sqr(g1) - 12
      *Lambda3*M112*Sqr(g1d) - 6*Lambda4*M112*Sqr(g1d) + 18*Lambda7*M122*Sqr(
      g1d) - 4*Lambda3*M112*Sqr(g1dp) - 2*Lambda4*M112*Sqr(g1dp) + 6*Lambda7*
      M122*Sqr(g1dp) + 24*Lambda3*M112*Sqr(g2) + 12*Lambda4*M112*Sqr(g2) - 72*
      Lambda7*M122*Sqr(g2) + 72*Lambda2*M222*Sqr(g2) + 11.25*M222*traceYuAdjYu*
      Sqr(g2) + 1.125*M222*Sqr(g1)*Sqr(g2) + 18*Lambda7*M122*Sqr(g2u) - 36*
      Lambda2*M222*Sqr(g2u) + 1.125*M222*Sqr(g1)*Sqr(g2u) - 2.25*M222*Sqr(g1d)*
      Sqr(g2u) + 20.625*M222*Sqr(g2)*Sqr(g2u) + 0.5*g2up*Conj(MassB)*(11*MassB*
      Cube(g2up) + 12*g1dp*Lambda7*Mu + 6*g2up*MassB*Sqr(g1dp) + 9*g2up*MassB*
      Sqr(g2u) + 6*g2up*MassWB*Sqr(g2u)) + 6*Lambda7*M122*Sqr(g2up) - 12*
      Lambda2*M222*Sqr(g2up) + 0.375*M222*Sqr(g1)*Sqr(g2up) - 0.75*M222*Sqr(
      g1dp)*Sqr(g2up) + 1.875*M222*Sqr(g2)*Sqr(g2up) - 2.25*M222*Sqr(g2u)*Sqr(
      g2up) + 1.5*Conj(MassWB)*(12*g1d*g2u*Lambda7*Mu - 24*MassWB*Quad(g2) + 13
      *MassWB*Quad(g2u) + Sqr(g2u)*(6*MassWB*(Sqr(g1d) - 4*Sqr(g2)) + (2*MassB
      + 3*MassWB)*Sqr(g2up))) + 40*M222*traceYuAdjYu*Sqr(g3) - 60*M222*Sqr(
      Lambda2) - 8*M112*Sqr(Lambda3) - 2*M222*Sqr(Lambda3) - 8*M112*Sqr(Lambda4
      ) - 2*M222*Sqr(Lambda4) - 12*M112*Sqr(Lambda5) - 3*M222*Sqr(Lambda5) - 18
      *M112*Sqr(Lambda6) + 3*M222*Sqr(Lambda6) - 18*M112*Sqr(Lambda7) - 27*M222
      *Sqr(Lambda7) + 6*g1d*g1dp*g2u*g2up*Sqr(Mu) - 2.16*Quad(g1)*Sqr(Mu) - 18*
      Quad(g2)*Sqr(Mu) + 18*Quad(g2u)*Sqr(Mu) + 4*Quad(g2up)*Sqr(Mu) - 1.8*Sqr(
      g1)*Sqr(g2u)*Sqr(Mu) + 10.5*Sqr(g1d)*Sqr(g2u)*Sqr(Mu) + 1.5*Sqr(g1dp)*Sqr
      (g2u)*Sqr(Mu) - 21*Sqr(g2)*Sqr(g2u)*Sqr(Mu) - 0.6*Sqr(g1)*Sqr(g2up)*Sqr(
      Mu) + 1.5*Sqr(g1d)*Sqr(g2up)*Sqr(Mu) + 4.5*Sqr(g1dp)*Sqr(g2up)*Sqr(Mu) -
      3*Sqr(g2)*Sqr(g2up)*Sqr(Mu) + 6*Sqr(g2u)*Sqr(g2up)*Sqr(Mu)));


   return beta_M222;
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


   return beta_M222;
}

} // namespace flexiblesusy
