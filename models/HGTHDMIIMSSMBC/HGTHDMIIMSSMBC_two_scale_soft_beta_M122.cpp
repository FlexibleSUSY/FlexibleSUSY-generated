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

// File generated at Sat 15 Oct 2016 15:25:18

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
 * Calculates the one-loop beta function of M122.
 *
 * @return one-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M122_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;


   double beta_M122;

   beta_M122 = Re(oneOver16PiSqr*(-6*Lambda6*M112 + 2*Lambda3*M122 + 4*
      Lambda4*M122 + 6*Lambda5*M122 - 6*Lambda7*M222 + 3*M122*traceYdAdjYd +
      M122*traceYeAdjYe + 3*M122*traceYuAdjYu + g1dp*g2up*MassB*Mu + 3*g1d*g2u*
      MassWB*Mu + g1dp*g2up*Conj(MassB)*Mu + 3*g1d*g2u*Conj(MassWB)*Mu - 0.9*
      M122*Sqr(g1) + 1.5*M122*Sqr(g1d) + 0.5*M122*Sqr(g1dp) - 4.5*M122*Sqr(g2)
      + 1.5*M122*Sqr(g2u) + 0.5*M122*Sqr(g2up)));


   return beta_M122;
}

/**
 * Calculates the two-loop beta function of M122.
 *
 * @return two-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M122_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;


   double beta_M122;

   const double beta_M122_1 = Re(0.0025*twoLoop*(1497*Power(g1,4)*M122 +
      5*Sqr(g1)*(-576*Lambda6*M112 + 192*Lambda3*M122 + 384*Lambda4*M122 + 576*
      Lambda5*M122 - 576*Lambda7*M222 + 50*M122*traceYdAdjYd + 150*M122*
      traceYeAdjYe + 170*M122*traceYuAdjYu + 45*M122*Sqr(g1d) + 15*M122*Sqr(
      g1dp) + 90*M122*Sqr(g2) + 45*M122*Sqr(g2u) + 15*M122*Sqr(g2up)) - 25*(
      -576*Lambda1*Lambda6*M112 - 192*Lambda3*Lambda6*M112 - 192*Lambda4*
      Lambda6*M112 - 192*Lambda5*Lambda6*M112 - 96*Lambda3*Lambda7*M112 - 96*
      Lambda4*Lambda7*M112 - 96*Lambda5*Lambda7*M112 + 45*Power(g1d,4)*M122 + 9
      *Power(g1dp,4)*M122 + 123*Power(g2,4)*M122 + 45*Power(g2u,4)*M122 + 9*
      Power(g2up,4)*M122 + 192*Lambda1*Lambda3*M122 + 192*Lambda2*Lambda3*M122
      + 192*Lambda1*Lambda4*M122 + 192*Lambda2*Lambda4*M122 + 96*Lambda3*
      Lambda4*M122 + 192*Lambda1*Lambda5*M122 + 192*Lambda2*Lambda5*M122 + 192*
      Lambda3*Lambda5*M122 + 192*Lambda4*Lambda5*M122 + 576*Lambda6*Lambda7*
      M122 - 96*Lambda3*Lambda6*M222 - 96*Lambda4*Lambda6*M222 - 96*Lambda5*
      Lambda6*M222 - 576*Lambda2*Lambda7*M222 - 192*Lambda3*Lambda7*M222 - 192*
      Lambda4*Lambda7*M222 - 192*Lambda5*Lambda7*M222 - 576*Lambda6*M112*
      traceYdAdjYd + 96*Lambda3*M122*traceYdAdjYd + 192*Lambda4*M122*
      traceYdAdjYd + 288*Lambda5*M122*traceYdAdjYd + 108*M122*
      traceYdAdjYdYdAdjYd + 264*M122*traceYdAdjYuYuAdjYd - 192*Lambda6*M112*
      traceYeAdjYe + 32*Lambda3*M122*traceYeAdjYe + 64*Lambda4*M122*
      traceYeAdjYe + 96*Lambda5*M122*traceYeAdjYe + 36*M122*traceYeAdjYeYeAdjYe
      + 96*Lambda3*M122*traceYuAdjYu + 192*Lambda4*M122*traceYuAdjYu + 288*
      Lambda5*M122*traceYuAdjYu - 576*Lambda7*M222*traceYuAdjYu + 108*M122*
      traceYuAdjYuYuAdjYu + 24*g1d*g1dp*g2u*(-4*g2up*M122 + g1dp*MassB*Mu) +
      576*Lambda6*M112*Sqr(g2) - 192*Lambda3*M122*Sqr(g2) - 384*Lambda4*M122*
      Sqr(g2) - 576*Lambda5*M122*Sqr(g2) + 576*Lambda7*M222*Sqr(g2) - 90*M122*
      traceYdAdjYd*Sqr(g2) - 30*M122*traceYeAdjYe*Sqr(g2) - 90*M122*
      traceYuAdjYu*Sqr(g2) + 48*Lambda3*M122*Sqr(g2u) + 96*Lambda4*M122*Sqr(g2u
      ) + 144*Lambda5*M122*Sqr(g2u) - 288*Lambda7*M222*Sqr(g2u) - 165*M122*Sqr(
      g2)*Sqr(g2u) + 3*Sqr(g1d)*(-96*Lambda6*M112 + M122*(16*Lambda3 + 32*
      Lambda4 + 48*Lambda5 + 6*Sqr(g1dp) - 55*Sqr(g2) + 28*Sqr(g2u))) + Sqr(
      g1dp)*(-96*Lambda6*M112 + M122*(16*(Lambda3 + 2*Lambda4 + 3*Lambda5) - 15
      *Sqr(g2) - 4*Sqr(g2up))) + 16*Lambda3*M122*Sqr(g2up) + 32*Lambda4*M122*
      Sqr(g2up) + 48*Lambda5*M122*Sqr(g2up) - 96*Lambda7*M222*Sqr(g2up) - 15*
      M122*Sqr(g2)*Sqr(g2up) + 18*M122*Sqr(g2u)*Sqr(g2up) - 320*M122*
      traceYdAdjYd*Sqr(g3) - 320*M122*traceYuAdjYu*Sqr(g3) - 96*M122*Sqr(
      Lambda1) - 96*M122*Sqr(Lambda2) - 48*M122*Sqr(Lambda5) + 192*M122*Sqr(
      Lambda6) + 192*M122*Sqr(Lambda7))));
   const double beta_M122_2 = Re(0.05*twoLoop*Mu*(-40*Power(g1dp,3)*g2up*
      MassB - 30*g1dp*Power(g2up,3)*MassB - 20*g1dp*g2up*Lambda3*MassB - 40*
      g1dp*g2up*Lambda4*MassB - 60*g1dp*g2up*Lambda5*MassB - 120*Power(g1d,3)*
      g2u*MassWB - 30*g1d*Power(g2u,3)*MassWB - 60*g1d*g2u*Lambda3*MassWB - 120
      *g1d*g2u*Lambda4*MassWB - 180*g1d*g2u*Lambda5*MassWB + 3*g1dp*g2up*MassB*
      Sqr(g1) + 9*g1d*g2u*MassWB*Sqr(g1) - 30*g1dp*g2up*MassB*Sqr(g1d) - 30*
      g1dp*g2up*MassWB*Sqr(g1d) - 30*g1d*g2u*MassWB*Sqr(g1dp) + 15*g1dp*g2up*
      MassB*Sqr(g2) + 285*g1d*g2u*MassWB*Sqr(g2) - 30*g1dp*g2up*MassWB*Sqr(g2u)
      - 30*g1d*g2u*MassB*Sqr(g2up) + 3*Conj(MassWB)*(3*g1d*g2u*Sqr(g1) - 5*(2*
      Power(g1d,3)*g2u + 2*g1dp*g2up*Sqr(g1d) + 2*g1dp*g2up*Sqr(g2u) + g1d*g2u*
      (4*Lambda3 + 8*Lambda4 + 12*Lambda5 - 19*Sqr(g2) + 8*Sqr(g2u) + 2*Sqr(
      g2up)))) - Conj(MassB)*(30*g1d*g2u*(Sqr(g1dp) + Sqr(g2up)) + g1dp*g2up*(
      -3*Sqr(g1) + 5*(4*Lambda3 + 8*Lambda4 + 12*Lambda5 + 6*Sqr(g1dp) - 3*Sqr(
      g2) + 6*Sqr(g2u) + 8*Sqr(g2up))))));

   beta_M122 = beta_M122_1 + beta_M122_2;


   return beta_M122;
}

/**
 * Calculates the three-loop beta function of M122.
 *
 * @return three-loop beta function
 */
double HGTHDMIIMSSMBC_soft_parameters::calc_beta_M122_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_M122;

   beta_M122 = 0;


   return beta_M122;
}

} // namespace flexiblesusy
