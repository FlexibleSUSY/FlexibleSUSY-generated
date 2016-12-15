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

// File generated at Thu 15 Dec 2016 12:43:16

#include "E6SSMtower_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mHd2.
 *
 * @return one-loop beta function
 */
double E6SSMtower_soft_parameters::calc_beta_mHd2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_mHd2;

   beta_mHd2 = Re(oneOver16PiSqr*(-0.7745966692414834*g1*Tr11 -
      0.9486832980505138*gN*Tr14 + 6*traceconjTYdTpTYd + 2*traceconjTYeTpTYe +
      6*tracemd2YdAdjYd + 2*traceme2YeAdjYe + 2*traceml2AdjYeYe + 6*
      tracemq2AdjYdYd + 6*mHd2*traceYdAdjYd + 2*mHd2*traceYeAdjYe + 2*mHd2*
      AbsSqr(Lambdax) + 2*mHu2*AbsSqr(Lambdax) + 2*ms2*AbsSqr(Lambdax) + 2*
      AbsSqr(TLambdax) - 1.2*AbsSqr(MassB)*Sqr(g1) - 6*AbsSqr(MassWB)*Sqr(g2) -
      1.8*AbsSqr(MassBp)*Sqr(gN)));


   return beta_mHd2;
}

/**
 * Calculates the two-loop beta function of mHd2.
 *
 * @return two-loop beta function
 */
double E6SSMtower_soft_parameters::calc_beta_mHd2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd =
      TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd =
      TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd =
      TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd =
      TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe =
      TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe =
      TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjTYu =
      TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu =
      TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double tracemd2YdAdjYdYdAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd =
      TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe =
      TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe =
      TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu =
      TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd =
      TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu =
      TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHd2;

   const double beta_mHd2_1 = Re(0.01*twoLoop*(2*Conj(MassB)*Sqr(g1)*(40*
      traceAdjYdTYd - 120*traceAdjYeTYe - 80*MassB*traceYdAdjYd + 240*MassB*
      traceYeAdjYe + 1782*MassB*Sqr(g1) + 90*(2*MassB + MassWB)*Sqr(g2) - 18*
      MassB*Sqr(gN) - 9*MassBp*Sqr(gN)) + Conj(MassBp)*Sqr(gN)*(120*
      traceAdjYdTYd + 40*traceAdjYeTYe - 240*MassBp*traceYdAdjYd - 80*MassBp*
      traceYeAdjYe + 400*MassBp*AbsSqr(Lambdax) - 18*(MassB + 2*MassBp)*Sqr(g1)
      + 540*MassBp*Sqr(g2) + 270*MassWB*Sqr(g2) + 5319*MassBp*Sqr(gN)) + 10*(
      -320*(traceAdjYdTYd - 2*MassG*traceYdAdjYd)*Conj(MassG)*Sqr(g3) + 3*Conj(
      MassWB)*Sqr(g2)*(6*(MassB + 2*MassWB)*Sqr(g1) + 290*MassWB*Sqr(g2) + 9*(
      MassBp + 2*MassWB)*Sqr(gN)) + 2*(30*Power(g2,4)*Tr22 + 7.348469228349534*
      g1*gN*Tr2U114 + 7.348469228349534*g1*gN*Tr2U141 - 15.491933384829668*g1*
      Tr31 - 18.973665961010276*gN*Tr34 - 180*tracemd2YdAdjYdYdAdjYd - 30*
      tracemd2YdAdjYuYuAdjYd - 60*traceme2YeAdjYeYeAdjYe - 60*
      traceml2AdjYeYeAdjYeYe - 180*tracemq2AdjYdYdAdjYdYd - 30*
      tracemq2AdjYdYdAdjYuYu - 30*tracemq2AdjYuYuAdjYdYd - 30*
      tracemu2YuAdjYdYdAdjYu - 180*traceYdAdjTYdTYdAdjYd - 30*
      traceYdAdjTYuTYuAdjYd - 180*traceYdAdjYdTYdAdjTYd - 180*mHd2*
      traceYdAdjYdYdAdjYd - 30*traceYdAdjYuTYuAdjTYd - 30*mHd2*
      traceYdAdjYuYuAdjYd - 30*mHu2*traceYdAdjYuYuAdjYd - 60*
      traceYeAdjTYeTYeAdjYe - 60*traceYeAdjYeTYeAdjTYe - 60*mHd2*
      traceYeAdjYeYeAdjYe - 30*traceYuAdjTYdTYdAdjYu - 30*traceYuAdjYdTYdAdjTYu
      + 6*Tr2U111*Sqr(g1) - 4*traceconjTYdTpTYd*Sqr(g1) + 4*MassB*
      traceconjTYdTpYd*Sqr(g1) + 12*traceconjTYeTpTYe*Sqr(g1) - 12*MassB*
      traceconjTYeTpYe*Sqr(g1) - 4*tracemd2YdAdjYd*Sqr(g1) + 12*traceme2YeAdjYe
      *Sqr(g1) + 12*traceml2AdjYeYe*Sqr(g1) - 4*tracemq2AdjYdYd*Sqr(g1) - 4*
      mHd2*traceYdAdjYd*Sqr(g1) + 12*mHd2*traceYeAdjYe*Sqr(g1) + 160*
      traceconjTYdTpTYd*Sqr(g3) - 160*MassG*traceconjTYdTpYd*Sqr(g3) + 160*
      tracemd2YdAdjYd*Sqr(g3) + 160*tracemq2AdjYdYd*Sqr(g3) + 160*mHd2*
      traceYdAdjYd*Sqr(g3) + 9*Tr2U144*Sqr(gN) - 6*traceconjTYdTpTYd*Sqr(gN) +
      6*MassBp*traceconjTYdTpYd*Sqr(gN) - 2*traceconjTYeTpTYe*Sqr(gN) + 2*
      MassBp*traceconjTYeTpYe*Sqr(gN) - 6*tracemd2YdAdjYd*Sqr(gN) - 2*
      traceme2YeAdjYe*Sqr(gN) - 2*traceml2AdjYeYe*Sqr(gN) - 6*tracemq2AdjYdYd*
      Sqr(gN) - 6*mHd2*traceYdAdjYd*Sqr(gN) - 2*mHd2*traceYeAdjYe*Sqr(gN) + 10*
      AbsSqr(Lambdax)*(-3*traceconjTKappaTpTKappa - 2*
      traceconjTLambda12TpTLambda12 - 3*traceconjTYuTpTYu - 3*mHd2*
      traceKappaAdjKappa - 3*mHu2*traceKappaAdjKappa - 6*ms2*traceKappaAdjKappa
      - 3*traceKappaAdjKappaconjmDx2 - 3*traceKappaconjmDxbar2AdjKappa - 2*
      mHd2*traceLambda12AdjLambda12 - 2*mHu2*traceLambda12AdjLambda12 - 4*ms2*
      traceLambda12AdjLambda12 - 2*traceLambda12AdjLambda12conjmH2I2 - 2*
      tracemH1I2AdjLambda12Lambda12 - 3*tracemq2AdjYuYu - 3*tracemu2YuAdjYu - 3
      *mHd2*traceYuAdjYu - 6*mHu2*traceYuAdjYu - 3*ms2*traceYuAdjYu + (mHd2 +
      mHu2 + ms2)*Sqr(gN)) - 60*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax)))));
   const double beta_mHd2_2 = Re(-2*twoLoop*(6*(mHu2 + ms2)*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) + Conj(Lambdax)*(3*traceconjTKappaTpKappa + 2*
      traceconjTLambda12TpLambda12 + 3*traceconjTYuTpYu + 12*Conj(TLambdax)*
      Lambdax + Conj(MassBp)*Sqr(gN))*TLambdax + Conj(TLambdax)*(Lambdax*(3*
      traceAdjKappaTKappa + 2*traceAdjLambda12TLambda12 + 3*traceAdjYuTYu +
      MassBp*Sqr(gN)) + (3*traceKappaAdjKappa + 2*traceLambda12AdjLambda12 + 3*
      traceYuAdjYu - Sqr(gN))*TLambdax)));

   beta_mHd2 = beta_mHd2_1 + beta_mHd2_2;


   return beta_mHd2;
}

/**
 * Calculates the three-loop beta function of mHd2.
 *
 * @return three-loop beta function
 */
double E6SSMtower_soft_parameters::calc_beta_mHd2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

} // namespace flexiblesusy
