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

// File generated at Fri 8 Jan 2016 12:30:53

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
 * Calculates the one-loop beta function of ms2.
 *
 * @return one-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_ms2_one_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2;
   const double traceKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappa;
   const double traceLambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2;
   const double Tr14 = TRACE_STRUCT.Tr14;


   double beta_ms2;

   beta_ms2 = Re(oneOver16PiSqr*(1.5811388300841898*gN*Tr14 + 6*
      traceconjTKappaTpTKappa + 4*traceconjTLambda12TpTLambda12 + 6*ms2*
      traceKappaAdjKappa + 6*traceKappaAdjKappaconjmDx2 + 6*
      traceKappaconjmDxbar2AdjKappa + 4*ms2*traceLambda12AdjLambda12 + 4*
      traceLambda12AdjLambda12conjmH2I2 + 4*tracemH1I2AdjLambda12Lambda12 + 4*(
      mHd2 + mHu2 + ms2)*AbsSqr(Lambdax) + 4*AbsSqr(TLambdax) - 5*AbsSqr(MassBp
      )*Sqr(gN)));


   return beta_ms2;
}

/**
 * Calculates the two-loop beta function of ms2.
 *
 * @return two-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_ms2_two_loop(const Soft_traces& soft_traces) const
{
   const double traceKappaAdjKappa = TRACE_STRUCT.traceKappaAdjKappa;
   const double traceLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTKappaTpTKappa =
      TRACE_STRUCT.traceconjTKappaTpTKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceconjTLambda12TpTLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpTLambda12;
   const double tracemH1I2AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
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
   const double traceKappaAdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappa;
   const double traceKappaAdjKappaTKappaAdjTKappa =
      TRACE_STRUCT.traceKappaAdjKappaTKappaAdjTKappa;
   const double traceKappaAdjTKappaTKappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjTKappaTKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12;
   const double traceLambda12AdjLambda12TLambda12AdjTLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12TLambda12AdjTLambda12;
   const double traceLambda12AdjTLambda12TLambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjTLambda12TLambda12AdjLambda12;
   const double tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 =
      TRACE_STRUCT.tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12;
   const double traceKappaAdjKappaKappaAdjKappaconjmDx2 =
      TRACE_STRUCT.traceKappaAdjKappaKappaAdjKappaconjmDx2;
   const double traceKappaAdjKappaKappaconjmDxbar2AdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaKappaconjmDxbar2AdjKappa;
   const double traceKappaAdjKappaconjmDx2KappaAdjKappa =
      TRACE_STRUCT.traceKappaAdjKappaconjmDx2KappaAdjKappa;
   const double traceKappaconjmDxbar2AdjKappaKappaAdjKappa =
      TRACE_STRUCT.traceKappaconjmDxbar2AdjKappaKappaAdjKappa;
   const double traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 =
      TRACE_STRUCT.traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2;
   const double traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 =
      TRACE_STRUCT.traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   double beta_ms2;

   const double beta_ms2_1 = Re(0.05*twoLoop*(-16*(2*traceAdjKappaTKappa
      + 3*traceAdjLambda12TLambda12 - 4*MassB*traceKappaAdjKappa - 6*MassB*
      traceLambda12AdjLambda12 - 6*MassB*AbsSqr(Lambdax))*Conj(MassB)*Sqr(g1) +
      3*Conj(MassBp)*Sqr(gN)*(8*(3*traceAdjKappaTKappa + 2*
      traceAdjLambda12TLambda12 - 6*MassBp*traceKappaAdjKappa - 4*MassBp*
      traceLambda12AdjLambda12) - 32*MassBp*AbsSqr(Lambdax) + 1065*MassBp*Sqr(
      gN)) + 4*(31.622776601683796*gN*Tr34 - 60*
      traceKappaAdjKappaconjmDx2KappaAdjKappa - 120*ms2*
      traceKappaAdjKappaKappaAdjKappa - 60*
      traceKappaAdjKappaKappaAdjKappaconjmDx2 - 60*
      traceKappaAdjKappaKappaconjmDxbar2AdjKappa - 120*
      traceKappaAdjKappaTKappaAdjTKappa - 120*traceKappaAdjTKappaTKappaAdjKappa
      - 60*traceKappaconjmDxbar2AdjKappaKappaAdjKappa - 40*
      traceLambda12AdjLambda12conjmH2I2Lambda12AdjLambda12 - 80*ms2*
      traceLambda12AdjLambda12Lambda12AdjLambda12 - 40*
      traceLambda12AdjLambda12Lambda12AdjLambda12conjmH2I2 - 80*
      traceLambda12AdjLambda12TLambda12AdjTLambda12 - 80*
      traceLambda12AdjTLambda12TLambda12AdjLambda12 - 80*
      tracemH1I2AdjLambda12Lambda12AdjLambda12Lambda12 - 60*traceconjTYdTpTYd*
      AbsSqr(Lambdax) - 20*traceconjTYeTpTYe*AbsSqr(Lambdax) - 60*
      traceconjTYuTpTYu*AbsSqr(Lambdax) - 60*tracemd2YdAdjYd*AbsSqr(Lambdax) -
      20*traceme2YeAdjYe*AbsSqr(Lambdax) - 20*traceml2AdjYeYe*AbsSqr(Lambdax) -
      60*tracemq2AdjYdYd*AbsSqr(Lambdax) - 60*tracemq2AdjYuYu*AbsSqr(Lambdax)
      - 60*tracemu2YuAdjYu*AbsSqr(Lambdax) - 120*mHd2*traceYdAdjYd*AbsSqr(
      Lambdax) - 60*mHu2*traceYdAdjYd*AbsSqr(Lambdax) - 60*ms2*traceYdAdjYd*
      AbsSqr(Lambdax) - 40*mHd2*traceYeAdjYe*AbsSqr(Lambdax) - 20*mHu2*
      traceYeAdjYe*AbsSqr(Lambdax) - 20*ms2*traceYeAdjYe*AbsSqr(Lambdax) - 60*
      mHd2*traceYuAdjYu*AbsSqr(Lambdax) - 120*mHu2*traceYuAdjYu*AbsSqr(Lambdax)
      - 60*ms2*traceYuAdjYu*AbsSqr(Lambdax) - 60*traceAdjYdTYd*Conj(TLambdax)*
      Lambdax - 20*traceAdjYeTYe*Conj(TLambdax)*Lambdax - 60*traceAdjYuTYu*Conj
      (TLambdax)*Lambdax - 8*MassB*traceconjTKappaTpKappa*Sqr(g1) + 8*
      traceconjTKappaTpTKappa*Sqr(g1) - 12*MassB*traceconjTLambda12TpLambda12*
      Sqr(g1) + 12*traceconjTLambda12TpTLambda12*Sqr(g1) + 8*ms2*
      traceKappaAdjKappa*Sqr(g1) + 8*traceKappaAdjKappaconjmDx2*Sqr(g1) + 8*
      traceKappaconjmDxbar2AdjKappa*Sqr(g1) + 12*ms2*traceLambda12AdjLambda12*
      Sqr(g1) + 12*traceLambda12AdjLambda12conjmH2I2*Sqr(g1) + 12*
      tracemH1I2AdjLambda12Lambda12*Sqr(g1) + 12*mHd2*AbsSqr(Lambdax)*Sqr(g1) +
      12*mHu2*AbsSqr(Lambdax)*Sqr(g1) + 12*ms2*AbsSqr(Lambdax)*Sqr(g1) - 12*
      MassB*Conj(TLambdax)*Lambdax*Sqr(g1) - 60*MassWB*
      traceconjTLambda12TpLambda12*Sqr(g2) + 60*traceconjTLambda12TpTLambda12*
      Sqr(g2) + 60*ms2*traceLambda12AdjLambda12*Sqr(g2) + 60*
      traceLambda12AdjLambda12conjmH2I2*Sqr(g2) + 60*
      tracemH1I2AdjLambda12Lambda12*Sqr(g2) + 60*mHd2*AbsSqr(Lambdax)*Sqr(g2) +
      60*mHu2*AbsSqr(Lambdax)*Sqr(g2) + 60*ms2*AbsSqr(Lambdax)*Sqr(g2) - 60*(
      traceAdjLambda12TLambda12 - 2*MassWB*traceLambda12AdjLambda12 - 2*MassWB*
      AbsSqr(Lambdax))*Conj(MassWB)*Sqr(g2) - 60*MassWB*Conj(TLambdax)*Lambdax*
      Sqr(g2) - 160*MassG*traceconjTKappaTpKappa*Sqr(g3) + 160*
      traceconjTKappaTpTKappa*Sqr(g3) + 160*ms2*traceKappaAdjKappa*Sqr(g3) +
      160*traceKappaAdjKappaconjmDx2*Sqr(g3) + 160*
      traceKappaconjmDxbar2AdjKappa*Sqr(g3) - 160*(traceAdjKappaTKappa - 2*
      MassG*traceKappaAdjKappa)*Conj(MassG)*Sqr(g3) + 25*Tr2U144*Sqr(gN) + 18*
      MassBp*traceconjTKappaTpKappa*Sqr(gN) - 18*traceconjTKappaTpTKappa*Sqr(gN
      ) + 12*MassBp*traceconjTLambda12TpLambda12*Sqr(gN) - 12*
      traceconjTLambda12TpTLambda12*Sqr(gN) - 18*ms2*traceKappaAdjKappa*Sqr(gN)
      - 18*traceKappaAdjKappaconjmDx2*Sqr(gN) - 18*
      traceKappaconjmDxbar2AdjKappa*Sqr(gN) - 12*ms2*traceLambda12AdjLambda12*
      Sqr(gN) - 12*traceLambda12AdjLambda12conjmH2I2*Sqr(gN) - 12*
      tracemH1I2AdjLambda12Lambda12*Sqr(gN) - 12*mHd2*AbsSqr(Lambdax)*Sqr(gN) -
      12*mHu2*AbsSqr(Lambdax)*Sqr(gN) - 12*ms2*AbsSqr(Lambdax)*Sqr(gN) + 12*
      MassBp*Conj(TLambdax)*Lambdax*Sqr(gN) - 80*mHd2*Sqr(Conj(Lambdax))*Sqr(
      Lambdax) - 80*mHu2*Sqr(Conj(Lambdax))*Sqr(Lambdax) - 80*ms2*Sqr(Conj(
      Lambdax))*Sqr(Lambdax) - 60*traceconjTYdTpYd*Conj(Lambdax)*TLambdax - 20*
      traceconjTYeTpYe*Conj(Lambdax)*TLambdax - 60*traceconjTYuTpYu*Conj(
      Lambdax)*TLambdax)));
   const double beta_ms2_2 = Re(-0.8*twoLoop*(15*traceYdAdjYd*Conj(
      TLambdax) + 5*traceYeAdjYe*Conj(TLambdax) + 15*traceYuAdjYu*Conj(TLambdax
      ) + 40*AbsSqr(Lambdax)*Conj(TLambdax) + 3*Conj(MassB)*Conj(Lambdax)*Sqr(
      g1) - 3*Conj(TLambdax)*Sqr(g1) + 15*Conj(MassWB)*Conj(Lambdax)*Sqr(g2) -
      15*Conj(TLambdax)*Sqr(g2) - 3*Conj(MassBp)*Conj(Lambdax)*Sqr(gN) + 3*Conj
      (TLambdax)*Sqr(gN))*TLambdax);

   beta_ms2 = beta_ms2_1 + beta_ms2_2;


   return beta_ms2;
}

/**
 * Calculates the three-loop beta function of ms2.
 *
 * @return three-loop beta function
 */
double E6SSM_soft_parameters::calc_beta_ms2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   double beta_ms2;

   beta_ms2 = 0;


   return beta_ms2;
}

} // namespace flexiblesusy
