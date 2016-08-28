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

// File generated at Sun 28 Aug 2016 15:10:11

#include "E6SSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mHu2.
 *
 * @return one-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_mHu2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHu2;

   beta_mHu2 = Re(oneOver16PiSqr*(0.7745966692414834*g1*Tr11 -
      0.6324555320336759*gN*Tr14 + 6*traceconjTYuTpTYu + 6*tracemq2AdjYuYu + 6*
      tracemu2YuAdjYu + 6*mHu2*traceYuAdjYu + 2*mHd2*AbsSqr(Lambdax) + 2*mHu2*
      AbsSqr(Lambdax) + 2*ms2*AbsSqr(Lambdax) + 2*AbsSqr(TLambdax) - 1.2*AbsSqr
      (MassB)*Sqr(g1) - 6*AbsSqr(MassWB)*Sqr(g2) - 0.8*AbsSqr(MassBp)*Sqr(gN)))
      ;


   return beta_mHu2;
}

/**
 * Calculates the two-loop beta function of mHu2.
 *
 * @return two-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_mHu2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjYuYuAdjYu = TRACE_STRUCT.traceYuAdjYuYuAdjYu;
   const double traceYuAdjYuTYuAdjTYu =
      TRACE_STRUCT.traceYuAdjYuTYuAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double traceYuAdjTYuTYuAdjYu =
      TRACE_STRUCT.traceYuAdjTYuTYuAdjYu;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemq2AdjYuYuAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYuYu;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double tracemu2YuAdjYuYuAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYuYuAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHu2;

   const double beta_mHu2_1 = Re(0.04*twoLoop*(Conj(MassB)*Sqr(g1)*(-40*
      traceAdjYuTYu + 80*MassB*traceYuAdjYu + 891*MassB*Sqr(g1) + 45*(2*MassB +
      MassWB)*Sqr(g2) + 36*MassB*Sqr(gN) + 18*MassBp*Sqr(gN)) + 3*Conj(MassBp)
      *Sqr(gN)*(5*traceAdjYuTYu - 10*MassBp*traceYuAdjYu + 6*(MassB + 2*MassBp)
      *Sqr(g1) + 20*MassBp*Sqr(g2) + 10*MassWB*Sqr(g2) + 192*MassBp*Sqr(gN) +
      25*Conj(Lambdax)*(2*MassBp*Lambdax - TLambdax)) + 5*(30*Power(g2,4)*Tr22
      - 4.898979485566356*g1*gN*Tr2U114 - 4.898979485566356*g1*gN*Tr2U141 +
      15.491933384829668*g1*Tr31 - 12.649110640673518*gN*Tr34 - 30*
      tracemd2YdAdjYuYuAdjYd - 30*tracemq2AdjYdYdAdjYuYu - 30*
      tracemq2AdjYuYuAdjYdYd - 180*tracemq2AdjYuYuAdjYuYu - 30*
      tracemu2YuAdjYdYdAdjYu - 180*tracemu2YuAdjYuYuAdjYu - 30*
      traceYdAdjTYuTYuAdjYd - 30*traceYdAdjYuTYuAdjTYd - 30*mHd2*
      traceYdAdjYuYuAdjYd - 30*mHu2*traceYdAdjYuYuAdjYd - 30*
      traceYuAdjTYdTYdAdjYu - 180*traceYuAdjTYuTYuAdjYu - 30*
      traceYuAdjYdTYdAdjTYu - 180*traceYuAdjYuTYuAdjTYu - 180*mHu2*
      traceYuAdjYuYuAdjYu - 30*traceconjTKappaTpTKappa*AbsSqr(Lambdax) - 20*
      traceconjTLambda12TpTLambda12*AbsSqr(Lambdax) - 30*traceconjTYdTpTYd*
      AbsSqr(Lambdax) - 10*traceconjTYeTpTYe*AbsSqr(Lambdax) - 30*mHd2*
      traceKappaAdjKappa*AbsSqr(Lambdax) - 30*mHu2*traceKappaAdjKappa*AbsSqr(
      Lambdax) - 60*ms2*traceKappaAdjKappa*AbsSqr(Lambdax) - 30*
      traceKappaAdjKappaconjmDx2*AbsSqr(Lambdax) - 30*
      traceKappaconjmDxbar2AdjKappa*AbsSqr(Lambdax) - 20*mHd2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 20*mHu2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 40*ms2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 20*
      traceLambda12AdjLambda12conjmH2I2*AbsSqr(Lambdax) - 30*tracemd2YdAdjYd*
      AbsSqr(Lambdax) - 10*traceme2YeAdjYe*AbsSqr(Lambdax) - 20*
      tracemH1I2AdjLambda12Lambda12*AbsSqr(Lambdax) - 10*traceml2AdjYeYe*AbsSqr
      (Lambdax) - 30*tracemq2AdjYdYd*AbsSqr(Lambdax) - 60*mHd2*traceYdAdjYd*
      AbsSqr(Lambdax) - 30*mHu2*traceYdAdjYd*AbsSqr(Lambdax) - 30*ms2*
      traceYdAdjYd*AbsSqr(Lambdax) - 20*mHd2*traceYeAdjYe*AbsSqr(Lambdax) - 10*
      mHu2*traceYeAdjYe*AbsSqr(Lambdax) - 10*ms2*traceYeAdjYe*AbsSqr(Lambdax) -
      30*traceAdjKappaTKappa*Conj(TLambdax)*Lambdax - 20*
      traceAdjLambda12TLambda12*Conj(TLambdax)*Lambdax - 30*traceAdjYdTYd*Conj(
      TLambdax)*Lambdax - 10*traceAdjYeTYe*Conj(TLambdax)*Lambdax + 6*Tr2U111*
      Sqr(g1) + 8*traceconjTYuTpTYu*Sqr(g1) - 8*MassB*traceconjTYuTpYu*Sqr(g1)
      + 8*tracemq2AdjYuYu*Sqr(g1) + 8*tracemu2YuAdjYu*Sqr(g1) + 8*mHu2*
      traceYuAdjYu*Sqr(g1) + 160*traceconjTYuTpTYu*Sqr(g3) - 160*MassG*
      traceconjTYuTpYu*Sqr(g3) + 160*tracemq2AdjYuYu*Sqr(g3) + 160*
      tracemu2YuAdjYu*Sqr(g3) + 160*mHu2*traceYuAdjYu*Sqr(g3) - 160*(
      traceAdjYuTYu - 2*MassG*traceYuAdjYu)*Conj(MassG)*Sqr(g3) + 4*Tr2U144*Sqr
      (gN) - 3*traceconjTYuTpTYu*Sqr(gN) + 3*MassBp*traceconjTYuTpYu*Sqr(gN) -
      3*tracemq2AdjYuYu*Sqr(gN) - 3*tracemu2YuAdjYu*Sqr(gN) - 3*mHu2*
      traceYuAdjYu*Sqr(gN) + 15*mHd2*AbsSqr(Lambdax)*Sqr(gN) + 15*mHu2*AbsSqr(
      Lambdax)*Sqr(gN) + 15*ms2*AbsSqr(Lambdax)*Sqr(gN) + 15*AbsSqr(TLambdax)*
      Sqr(gN) - 15*MassBp*Conj(TLambdax)*Lambdax*Sqr(gN) + 3*Conj(MassWB)*Sqr(
      g2)*(3*(MassB + 2*MassWB)*Sqr(g1) + 145*MassWB*Sqr(g2) + 2*(MassBp + 2*
      MassWB)*Sqr(gN)) - 60*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 60*mHu2*Sqr(
      Conj(Lambdax))*Sqr(Lambdax) - 60*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 30
      *traceconjTKappaTpKappa*Conj(Lambdax)*TLambdax - 20*
      traceconjTLambda12TpLambda12*Conj(Lambdax)*TLambdax - 30*traceconjTYdTpYd
      *Conj(Lambdax)*TLambdax - 10*traceconjTYeTpYe*Conj(Lambdax)*TLambdax)));
   const double beta_mHu2_2 = Re(-2*twoLoop*(3*traceKappaAdjKappa + 2*
      traceLambda12AdjLambda12 + 3*traceYdAdjYd + traceYeAdjYe + 12*AbsSqr(
      Lambdax))*AbsSqr(TLambdax));

   beta_mHu2 = beta_mHu2_1 + beta_mHu2_2;


   return beta_mHu2;
}

/**
 * Calculates the three-loop beta function of mHu2.
 *
 * @return three-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_mHu2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHu2;

   beta_mHu2 = 0;


   return beta_mHu2;
}

} // namespace flexiblesusy
