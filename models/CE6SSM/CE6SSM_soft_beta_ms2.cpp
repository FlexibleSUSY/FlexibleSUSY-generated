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

// File generated at Sun 4 Aug 2019 17:11:58

#include "CE6SSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of ms2.
 *
 * @return 1-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_ms2_1_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12;
   const double traceconjTKappaTpTKappa = TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpTLambda12 = TRACE_STRUCT.
      traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 = TRACE_STRUCT.
      tracemH1I2AdjLambda12Lambda12;
   const double traceKappaAdjKappaconjmDx2 = TRACE_STRUCT.
      traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa = TRACE_STRUCT.
      traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 = TRACE_STRUCT.
      traceLambda12AdjLambda12conjmH2I2;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_ms2;

   beta_ms2 = Re(0.5*oneOver16PiSqr*(3.1622776601683795*gN*Tr14 + 12*
      traceconjTKappaTpTKappa + 8*traceconjTLambda12TpTLambda12 + 12*ms2*
      traceKappaAdjKappa + 12*traceKappaAdjKappaconjmDx2 + 12*
      traceKappaconjmDxbar2AdjKappa + 8*ms2*traceLambda12AdjLambda12 + 8*
      traceLambda12AdjLambda12conjmH2I2 + 8*tracemH1I2AdjLambda12Lambda12 + 8*
      mHd2*AbsSqr(Lambdax) + 8*mHu2*AbsSqr(Lambdax) + 8*ms2*AbsSqr(Lambdax) + 8
      *AbsSqr(TLambdax) - 10*AbsSqr(MassBp)*Sqr(gN)));


   return beta_ms2;
}

/**
 * Calculates the 2-loop beta function of ms2.
 *
 * @return 2-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_ms2_2_loop(const Soft_traces& soft_traces) const
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
   const double traceKappaAdjKappaKappaAdjKappa = TRACE_STRUCT.
      traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjTKappa = TRACE_STRUCT.
      traceKappaAdjKappaTKappaAdjTKappa;
   const double traceKappaAdjTKappaTKappaAdjKappa = TRACE_STRUCT.
      traceKappaAdjTKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12TLambda12AdjTLambda12 = TRACE_STRUCT.
      traceLambda12AdjLambda12TLambda12AdjTLambda12;
   const double traceLambda12AdjTLambda12TLambda12AdjLambda12 = TRACE_STRUCT.
      traceLambda12AdjTLambda12TLambda12AdjLambda12;
   const double tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 = TRACE_STRUCT
      .tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12;
   const double traceKappaAdjKappaKappaAdjKappaconjmDx2 = TRACE_STRUCT.
      traceKappaAdjKappaKappaAdjKappaconjmDx2;
   const double traceKappaAdjKappaKappaconjmDxbar2AdjKappa = TRACE_STRUCT.
      traceKappaAdjKappaKappaconjmDxbar2AdjKappa;
   const double traceKappaAdjKappaconjmDx2KappaAdjKappa = TRACE_STRUCT.
      traceKappaAdjKappaconjmDx2KappaAdjKappa;
   const double traceKappaconjmDxbar2AdjKappaKappaAdjKappa = TRACE_STRUCT.
      traceKappaconjmDxbar2AdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2;
   const double traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_ms2;

   const double beta_ms2_1 = Re(0.05*twoLoop*(126.49110640673518*gN*Tr34 - 240*
      traceKappaAdjKappaconjmDx2KappaAdjKappa - 480*ms2*
      traceKappaAdjKappaKappaAdjKappa - 240*
      traceKappaAdjKappaKappaAdjKappaconjmDx2 - 240*
      traceKappaAdjKappaKappaconjmDxbar2AdjKappa - 480*
      traceKappaAdjKappaTKappaAdjTKappa - 480*traceKappaAdjTKappaTKappaAdjKappa
       - 240*traceKappaconjmDxbar2AdjKappaKappaAdjKappa - 160*
      traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 - 320*ms2*
      traceLambda12AdjLambda12Lambda12AdjLambda12 - 160*
      traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 - 320*
      traceLambda12AdjLambda12TLambda12AdjTLambda12 - 320*
      traceLambda12AdjTLambda12TLambda12AdjLambda12 - 320*
      tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 - 240*traceconjTYdTpTYd*
      AbsSqr(Lambdax) - 80*traceconjTYeTpTYe*AbsSqr(Lambdax) - 240*
      traceconjTYuTpTYu*AbsSqr(Lambdax) - 240*tracemd2YdAdjYd*AbsSqr(Lambdax) -
      80*traceme2YeAdjYe*AbsSqr(Lambdax) - 80*traceml2AdjYeYe*AbsSqr(Lambdax) -
      240*tracemq2AdjYdYd*AbsSqr(Lambdax) - 240*tracemq2AdjYuYu*AbsSqr(Lambdax)
      - 240*tracemu2YuAdjYu*AbsSqr(Lambdax) - 480*mHd2*traceYdAdjYd*AbsSqr(
      Lambdax) - 240*mHu2*traceYdAdjYd*AbsSqr(Lambdax) - 240*ms2*traceYdAdjYd*
      AbsSqr(Lambdax) - 160*mHd2*traceYeAdjYe*AbsSqr(Lambdax) - 80*mHu2*
      traceYeAdjYe*AbsSqr(Lambdax) - 80*ms2*traceYeAdjYe*AbsSqr(Lambdax) - 240*
      mHd2*traceYuAdjYu*AbsSqr(Lambdax) - 480*mHu2*traceYuAdjYu*AbsSqr(Lambdax)
      - 240*ms2*traceYuAdjYu*AbsSqr(Lambdax) - 240*traceAdjYdTYd*Conj(TLambdax)
      *Lambdax - 80*traceAdjYeTYe*Conj(TLambdax)*Lambdax - 240*traceAdjYuTYu*
      Conj(TLambdax)*Lambdax + 3195*AbsSqr(MassBp)*Quad(gN) - 32*MassB*
      traceconjTKappaTpKappa*Sqr(g1) + 32*traceconjTKappaTpTKappa*Sqr(g1) - 48*
      MassB*traceconjTLambda12TpLambda12*Sqr(g1) + 48*
      traceconjTLambda12TpTLambda12*Sqr(g1) + 32*ms2*traceKappaAdjKappa*Sqr(g1)
      + 32*traceKappaAdjKappaconjmDx2*Sqr(g1) + 32*
      traceKappaconjmDxbar2AdjKappa*Sqr(g1) + 48*ms2*traceLambda12AdjLambda12*
      Sqr(g1) + 48*traceLambda12AdjLambda12conjmH2I2*Sqr(g1) + 48*
      tracemH1I2AdjLambda12Lambda12*Sqr(g1) + 64*traceKappaAdjKappa*AbsSqr(
      MassB)*Sqr(g1) + 96*traceLambda12AdjLambda12*AbsSqr(MassB)*Sqr(g1) + 48*
      mHd2*AbsSqr(Lambdax)*Sqr(g1) + 48*mHu2*AbsSqr(Lambdax)*Sqr(g1) + 48*ms2*
      AbsSqr(Lambdax)*Sqr(g1) + 96*AbsSqr(MassB)*AbsSqr(Lambdax)*Sqr(g1) - 32*
      traceAdjKappaTKappa*Conj(MassB)*Sqr(g1) - 48*traceAdjLambda12TLambda12*
      Conj(MassB)*Sqr(g1) - 48*MassB*Conj(TLambdax)*Lambdax*Sqr(g1) - 240*
      MassWB*traceconjTLambda12TpLambda12*Sqr(g2) + 240*
      traceconjTLambda12TpTLambda12*Sqr(g2) + 240*ms2*traceLambda12AdjLambda12*
      Sqr(g2) + 240*traceLambda12AdjLambda12conjmH2I2*Sqr(g2) + 240*
      tracemH1I2AdjLambda12Lambda12*Sqr(g2) + 480*traceLambda12AdjLambda12*
      AbsSqr(MassWB)*Sqr(g2) + 240*mHd2*AbsSqr(Lambdax)*Sqr(g2) + 240*mHu2*
      AbsSqr(Lambdax)*Sqr(g2) + 240*ms2*AbsSqr(Lambdax)*Sqr(g2) + 480*AbsSqr(
      MassWB)*AbsSqr(Lambdax)*Sqr(g2) - 240*traceAdjLambda12TLambda12*Conj(
      MassWB)*Sqr(g2) - 240*MassWB*Conj(TLambdax)*Lambdax*Sqr(g2) - 640*MassG*
      traceconjTKappaTpKappa*Sqr(g3) + 640*traceconjTKappaTpTKappa*Sqr(g3) +
      640*ms2*traceKappaAdjKappa*Sqr(g3) + 640*traceKappaAdjKappaconjmDx2*Sqr(
      g3) + 640*traceKappaconjmDxbar2AdjKappa*Sqr(g3) + 1280*traceKappaAdjKappa
      *AbsSqr(MassG)*Sqr(g3) - 640*traceAdjKappaTKappa*Conj(MassG)*Sqr(g3) +
      100*Tr2U144*Sqr(gN) + 72*MassBp*traceconjTKappaTpKappa*Sqr(gN) - 72*
      traceconjTKappaTpTKappa*Sqr(gN) + 48*MassBp*traceconjTLambda12TpLambda12*
      Sqr(gN) - 48*traceconjTLambda12TpTLambda12*Sqr(gN) - 72*ms2*
      traceKappaAdjKappa*Sqr(gN) - 72*traceKappaAdjKappaconjmDx2*Sqr(gN) - 72*
      traceKappaconjmDxbar2AdjKappa*Sqr(gN) - 48*ms2*traceLambda12AdjLambda12*
      Sqr(gN) - 48*traceLambda12AdjLambda12conjmH2I2*Sqr(gN) - 48*
      tracemH1I2AdjLambda12Lambda12*Sqr(gN) - 144*traceKappaAdjKappa*AbsSqr(
      MassBp)*Sqr(gN) - 96*traceLambda12AdjLambda12*AbsSqr(MassBp)*Sqr(gN) - 48
      *mHd2*AbsSqr(Lambdax)*Sqr(gN) - 48*mHu2*AbsSqr(Lambdax)*Sqr(gN) - 48*ms2*
      AbsSqr(Lambdax)*Sqr(gN) - 96*AbsSqr(MassBp)*AbsSqr(Lambdax)*Sqr(gN) + 72*
      traceAdjKappaTKappa*Conj(MassBp)*Sqr(gN) + 48*traceAdjLambda12TLambda12*
      Conj(MassBp)*Sqr(gN) + 48*MassBp*Conj(TLambdax)*Lambdax*Sqr(gN) - 320*
      mHd2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 320*mHu2*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) - 320*ms2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 240*traceconjTYdTpYd
      *Conj(Lambdax)*TLambdax - 80*traceconjTYeTpYe*Conj(Lambdax)*TLambdax -
      240*traceconjTYuTpYu*Conj(Lambdax)*TLambdax));
   const double beta_ms2_2 = Re(-0.8*twoLoop*(15*traceYdAdjYd*Conj(TLambdax) +
      5*traceYeAdjYe*Conj(TLambdax) + 15*traceYuAdjYu*Conj(TLambdax) + 40*
      AbsSqr(Lambdax)*Conj(TLambdax) + 3*Conj(MassB)*Conj(Lambdax)*Sqr(g1) - 3*
      Conj(TLambdax)*Sqr(g1) + 15*Conj(MassWB)*Conj(Lambdax)*Sqr(g2) - 15*Conj(
      TLambdax)*Sqr(g2) - 3*Conj(MassBp)*Conj(Lambdax)*Sqr(gN) + 3*Conj(
      TLambdax)*Sqr(gN))*TLambdax);

   beta_ms2 = beta_ms2_1 + beta_ms2_2;


   return beta_ms2;
}

/**
 * Calculates the 3-loop beta function of ms2.
 *
 * @return 3-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_ms2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_ms2;

   beta_ms2 = 0;


   return beta_ms2;
}

/**
 * Calculates the 4-loop beta function of ms2.
 *
 * @return 4-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_ms2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_ms2;

   beta_ms2 = 0;


   return beta_ms2;
}

/**
 * Calculates the 5-loop beta function of ms2.
 *
 * @return 5-loop beta function
 */
double CE6SSM_soft_parameters::calc_beta_ms2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_ms2;

   beta_ms2 = 0;


   return beta_ms2;
}

} // namespace flexiblesusy
