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

// File generated at Fri 10 Apr 2020 20:02:12

#include "E6SSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mHd2.
 *
 * @return 1-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_mHd2_1_loop(const Soft_traces& soft_traces) const
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

   beta_mHd2 = Re(0.1*oneOver16PiSqr*(-7.745966692414834*g1*Tr11 -
      9.486832980505138*gN*Tr14 + 60*traceconjTYdTpTYd + 20*traceconjTYeTpTYe +
      60*tracemd2YdAdjYd + 20*traceme2YeAdjYe + 20*traceml2AdjYeYe + 60*
      tracemq2AdjYdYd + 60*mHd2*traceYdAdjYd + 20*mHd2*traceYeAdjYe + 20*mHd2*
      AbsSqr(Lambdax) + 20*mHu2*AbsSqr(Lambdax) + 20*ms2*AbsSqr(Lambdax) + 20*
      AbsSqr(TLambdax) - 12*AbsSqr(MassB)*Sqr(g1) - 60*AbsSqr(MassWB)*Sqr(g2) -
      18*AbsSqr(MassBp)*Sqr(gN)));


   return beta_mHd2;
}

/**
 * Calculates the 2-loop beta function of mHd2.
 *
 * @return 2-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_mHd2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 = TRACE_STRUCT.
      traceAdjLambda12TLambda12;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTKappaTpKappa = TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa = TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 = TRACE_STRUCT.
      traceconjTLambda12TpLambda12;
   const double traceconjTLambda12TpTLambda12 = TRACE_STRUCT.
      traceconjTLambda12TpTLambda12;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double tracemH1I2AdjLambda12Lambda12 = TRACE_STRUCT.
      tracemH1I2AdjLambda12Lambda12;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceKappaAdjKappaconjmDx2 = TRACE_STRUCT.
      traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa = TRACE_STRUCT.
      traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 = TRACE_STRUCT.
      traceLambda12AdjLambda12conjmH2I2;
   const double traceYdAdjYdYdAdjYd = TRACE_STRUCT.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYdTYdAdjTYd = TRACE_STRUCT.traceYdAdjYdTYdAdjTYd;
   const double traceYdAdjYuYuAdjYd = TRACE_STRUCT.traceYdAdjYuYuAdjYd;
   const double traceYdAdjYuTYuAdjTYd = TRACE_STRUCT.traceYdAdjYuTYuAdjTYd;
   const double traceYdAdjTYdTYdAdjYd = TRACE_STRUCT.traceYdAdjTYdTYdAdjYd;
   const double traceYdAdjTYuTYuAdjYd = TRACE_STRUCT.traceYdAdjTYuTYuAdjYd;
   const double traceYeAdjYeYeAdjYe = TRACE_STRUCT.traceYeAdjYeYeAdjYe;
   const double traceYeAdjYeTYeAdjTYe = TRACE_STRUCT.traceYeAdjYeTYeAdjTYe;
   const double traceYeAdjTYeTYeAdjYe = TRACE_STRUCT.traceYeAdjTYeTYeAdjYe;
   const double traceYuAdjYdTYdAdjTYu = TRACE_STRUCT.traceYuAdjYdTYdAdjTYu;
   const double traceYuAdjTYdTYdAdjYu = TRACE_STRUCT.traceYuAdjTYdTYdAdjYu;
   const double tracemd2YdAdjYdYdAdjYd = TRACE_STRUCT.tracemd2YdAdjYdYdAdjYd;
   const double tracemd2YdAdjYuYuAdjYd = TRACE_STRUCT.tracemd2YdAdjYuYuAdjYd;
   const double traceme2YeAdjYeYeAdjYe = TRACE_STRUCT.traceme2YeAdjYeYeAdjYe;
   const double traceml2AdjYeYeAdjYeYe = TRACE_STRUCT.traceml2AdjYeYeAdjYeYe;
   const double tracemq2AdjYdYdAdjYdYd = TRACE_STRUCT.tracemq2AdjYdYdAdjYdYd;
   const double tracemq2AdjYdYdAdjYuYu = TRACE_STRUCT.tracemq2AdjYdYdAdjYuYu;
   const double tracemq2AdjYuYuAdjYdYd = TRACE_STRUCT.tracemq2AdjYuYuAdjYdYd;
   const double tracemu2YuAdjYdYdAdjYu = TRACE_STRUCT.tracemu2YuAdjYdYdAdjYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_mHd2;

   const double beta_mHd2_1 = Re(0.01*twoLoop*(146.96938456699067*g1*gN*Tr2U114
       + 146.96938456699067*g1*gN*Tr2U141 - 309.83866769659335*g1*Tr31 -
      379.47331922020555*gN*Tr34 - 3600*tracemd2YdAdjYdYdAdjYd - 600*
      tracemd2YdAdjYuYuAdjYd - 1200*traceme2YeAdjYeYeAdjYe - 1200*
      traceml2AdjYeYeAdjYeYe - 3600*tracemq2AdjYdYdAdjYdYd - 600*
      tracemq2AdjYdYdAdjYuYu - 600*tracemq2AdjYuYuAdjYdYd - 600*
      tracemu2YuAdjYdYdAdjYu - 3600*traceYdAdjTYdTYdAdjYd - 600*
      traceYdAdjTYuTYuAdjYd - 3600*traceYdAdjYdTYdAdjTYd - 3600*mHd2*
      traceYdAdjYdYdAdjYd - 600*traceYdAdjYuTYuAdjTYd - 600*mHd2*
      traceYdAdjYuYuAdjYd - 600*mHu2*traceYdAdjYuYuAdjYd - 1200*
      traceYeAdjTYeTYeAdjYe - 1200*traceYeAdjYeTYeAdjTYe - 1200*mHd2*
      traceYeAdjYeYeAdjYe - 600*traceYuAdjTYdTYdAdjYu - 600*
      traceYuAdjYdTYdAdjTYu - 600*traceconjTKappaTpTKappa*AbsSqr(Lambdax) - 400
      *traceconjTLambda12TpTLambda12*AbsSqr(Lambdax) - 600*traceconjTYuTpTYu*
      AbsSqr(Lambdax) - 600*mHd2*traceKappaAdjKappa*AbsSqr(Lambdax) - 600*mHu2*
      traceKappaAdjKappa*AbsSqr(Lambdax) - 1200*ms2*traceKappaAdjKappa*AbsSqr(
      Lambdax) - 600*traceKappaAdjKappaconjmDx2*AbsSqr(Lambdax) - 600*
      traceKappaconjmDxbar2AdjKappa*AbsSqr(Lambdax) - 400*mHd2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 400*mHu2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 800*ms2*
      traceLambda12AdjLambda12*AbsSqr(Lambdax) - 400*
      traceLambda12AdjLambda12conjmH2I2*AbsSqr(Lambdax) - 400*
      tracemH1I2AdjLambda12Lambda12*AbsSqr(Lambdax) - 600*tracemq2AdjYuYu*
      AbsSqr(Lambdax) - 600*tracemu2YuAdjYu*AbsSqr(Lambdax) - 600*mHd2*
      traceYuAdjYu*AbsSqr(Lambdax) - 1200*mHu2*traceYuAdjYu*AbsSqr(Lambdax) -
      600*ms2*traceYuAdjYu*AbsSqr(Lambdax) + 3564*AbsSqr(MassB)*Quad(g1) + 600*
      Tr22*Quad(g2) + 8700*AbsSqr(MassWB)*Quad(g2) + 5319*AbsSqr(MassBp)*Quad(
      gN) + 120*Tr2U111*Sqr(g1) - 80*traceconjTYdTpTYd*Sqr(g1) + 80*MassB*
      traceconjTYdTpYd*Sqr(g1) + 240*traceconjTYeTpTYe*Sqr(g1) - 240*MassB*
      traceconjTYeTpYe*Sqr(g1) - 80*tracemd2YdAdjYd*Sqr(g1) + 240*
      traceme2YeAdjYe*Sqr(g1) + 240*traceml2AdjYeYe*Sqr(g1) - 80*
      tracemq2AdjYdYd*Sqr(g1) - 80*mHd2*traceYdAdjYd*Sqr(g1) + 240*mHd2*
      traceYeAdjYe*Sqr(g1) - 160*traceYdAdjYd*AbsSqr(MassB)*Sqr(g1) + 480*
      traceYeAdjYe*AbsSqr(MassB)*Sqr(g1) + 80*traceAdjYdTYd*Conj(MassB)*Sqr(g1)
      - 240*traceAdjYeTYe*Conj(MassB)*Sqr(g1) + 360*AbsSqr(MassB)*Sqr(g1)*Sqr(
      g2) + 360*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 180*MassWB*Conj(MassB)*Sqr(g1)
      *Sqr(g2) + 180*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 3200*
      traceconjTYdTpTYd*Sqr(g3) - 3200*MassG*traceconjTYdTpYd*Sqr(g3) + 3200*
      tracemd2YdAdjYd*Sqr(g3) + 3200*tracemq2AdjYdYd*Sqr(g3) + 3200*mHd2*
      traceYdAdjYd*Sqr(g3) + 6400*traceYdAdjYd*AbsSqr(MassG)*Sqr(g3) - 3200*
      traceAdjYdTYd*Conj(MassG)*Sqr(g3) + 180*Tr2U144*Sqr(gN) - 120*
      traceconjTYdTpTYd*Sqr(gN) + 120*MassBp*traceconjTYdTpYd*Sqr(gN) - 40*
      traceconjTYeTpTYe*Sqr(gN) + 40*MassBp*traceconjTYeTpYe*Sqr(gN) - 120*
      tracemd2YdAdjYd*Sqr(gN) - 40*traceme2YeAdjYe*Sqr(gN) - 40*traceml2AdjYeYe
      *Sqr(gN) - 120*tracemq2AdjYdYd*Sqr(gN) - 120*mHd2*traceYdAdjYd*Sqr(gN) -
      40*mHd2*traceYeAdjYe*Sqr(gN) - 240*traceYdAdjYd*AbsSqr(MassBp)*Sqr(gN) -
      80*traceYeAdjYe*AbsSqr(MassBp)*Sqr(gN) + 200*mHd2*AbsSqr(Lambdax)*Sqr(gN)
      + 200*mHu2*AbsSqr(Lambdax)*Sqr(gN) + 200*ms2*AbsSqr(Lambdax)*Sqr(gN) +
      400*AbsSqr(MassBp)*AbsSqr(Lambdax)*Sqr(gN) + 120*traceAdjYdTYd*Conj(
      MassBp)*Sqr(gN) + 40*traceAdjYeTYe*Conj(MassBp)*Sqr(gN) - 36*AbsSqr(MassB
      )*Sqr(g1)*Sqr(gN) - 36*AbsSqr(MassBp)*Sqr(g1)*Sqr(gN) - 18*MassBp*Conj(
      MassB)*Sqr(g1)*Sqr(gN) - 18*MassB*Conj(MassBp)*Sqr(g1)*Sqr(gN) + 540*
      AbsSqr(MassBp)*Sqr(g2)*Sqr(gN) + 540*AbsSqr(MassWB)*Sqr(g2)*Sqr(gN) + 270
      *MassWB*Conj(MassBp)*Sqr(g2)*Sqr(gN) + 270*MassBp*Conj(MassWB)*Sqr(g2)*
      Sqr(gN) - 1200*mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax)));
   const double beta_mHd2_2 = Re(-2*twoLoop*(3*traceKappaAdjKappa*AbsSqr(
      TLambdax) + 2*traceLambda12AdjLambda12*AbsSqr(TLambdax) + 3*traceYuAdjYu*
      AbsSqr(TLambdax) + 12*AbsSqr(Lambdax)*AbsSqr(TLambdax) + 3*
      traceAdjKappaTKappa*Conj(TLambdax)*Lambdax + 2*traceAdjLambda12TLambda12*
      Conj(TLambdax)*Lambdax + 3*traceAdjYuTYu*Conj(TLambdax)*Lambdax - AbsSqr(
      TLambdax)*Sqr(gN) + MassBp*Conj(TLambdax)*Lambdax*Sqr(gN) + 6*mHu2*Sqr(
      Conj(Lambdax))*Sqr(Lambdax) + 6*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax) + 3*
      traceconjTKappaTpKappa*Conj(Lambdax)*TLambdax + 2*
      traceconjTLambda12TpLambda12*Conj(Lambdax)*TLambdax + 3*traceconjTYuTpYu*
      Conj(Lambdax)*TLambdax + Conj(MassBp)*Conj(Lambdax)*Sqr(gN)*TLambdax));

   beta_mHd2 = beta_mHd2_1 + beta_mHd2_2;


   return beta_mHd2;
}

/**
 * Calculates the 3-loop beta function of mHd2.
 *
 * @return 3-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_mHd2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

/**
 * Calculates the 4-loop beta function of mHd2.
 *
 * @return 4-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_mHd2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

/**
 * Calculates the 5-loop beta function of mHd2.
 *
 * @return 5-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_mHd2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_mHd2;

   beta_mHd2 = 0;


   return beta_mHd2;
}

} // namespace flexiblesusy
