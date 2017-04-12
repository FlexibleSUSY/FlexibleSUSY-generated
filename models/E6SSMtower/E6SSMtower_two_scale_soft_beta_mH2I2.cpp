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

// File generated at Wed 12 Apr 2017 11:27:10

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
 * Calculates the one-loop beta function of mH2I2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMtower_soft_parameters::calc_beta_mH2I2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,2,2> beta_mH2I2;

   beta_mH2I2 = (oneOver16PiSqr*(2*ms2*(Lambda12.conjugate()*(Lambda12)
      .transpose()) + 2*(TLambda12.conjugate()*(TLambda12).transpose()) + mH2I2
      *Lambda12.conjugate()*(Lambda12).transpose() + 2*(Lambda12.conjugate()*
      mH1I2.conjugate()*(Lambda12).transpose()) + Lambda12.conjugate()*(
      Lambda12).transpose()*mH2I2 + 0.2*(3.872983346207417*g1*Tr11 -
      3.1622776601683795*gN*Tr14 - 6*AbsSqr(MassB)*Sqr(g1) - 30*AbsSqr(MassWB)*
      Sqr(g2) - 4*AbsSqr(MassBp)*Sqr(gN))*UNITMATRIX(2))).real();


   return beta_mH2I2;
}

/**
 * Calculates the two-loop beta function of mH2I2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMtower_soft_parameters::calc_beta_mH2I2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,2,2> beta_mH2I2;

   beta_mH2I2 = (twoLoop*((-6*traceconjTKappaTpTKappa - 4*
      traceconjTLambda12TpTLambda12 - 12*ms2*traceKappaAdjKappa - 6*
      traceKappaAdjKappaconjmDx2 - 6*traceKappaconjmDxbar2AdjKappa - 8*ms2*
      traceLambda12AdjLambda12 - 4*traceLambda12AdjLambda12conjmH2I2 - 4*
      tracemH1I2AdjLambda12Lambda12 - 4*(mHd2 + mHu2 + 2*ms2)*AbsSqr(Lambdax) -
      4*AbsSqr(TLambdax) + 3*ms2*Sqr(gN) + 6*AbsSqr(MassBp)*Sqr(gN))*(
      Lambda12.conjugate()*(Lambda12).transpose()) + (-6*traceconjTKappaTpKappa
      - 4*traceconjTLambda12TpLambda12 - 4*Conj(TLambdax)*Lambdax - 3*Conj(
      MassBp)*Sqr(gN))*(Lambda12.conjugate()*(TLambda12).transpose()) + (-6*
      traceAdjKappaTKappa - 4*traceAdjLambda12TLambda12 - 3*MassBp*Sqr(gN) - 4*
      Conj(Lambdax)*TLambdax)*(TLambda12.conjugate()*(Lambda12).transpose()) +
      (-6*traceKappaAdjKappa - 4*traceLambda12AdjLambda12 - 4*AbsSqr(Lambdax) +
      3*Sqr(gN))*(TLambda12.conjugate()*(TLambda12).transpose()) + (-3*
      traceKappaAdjKappa - 2*traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + 1.5
      *Sqr(gN))*(mH2I2*Lambda12.conjugate()*(Lambda12).transpose()) + (-6*
      traceKappaAdjKappa - 4*traceLambda12AdjLambda12 - 4*AbsSqr(Lambdax) + 3*
      Sqr(gN))*(Lambda12.conjugate()*mH1I2.conjugate()*(Lambda12).transpose())
      + (-3*traceKappaAdjKappa - 2*traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax)
      + 1.5*Sqr(gN))*(Lambda12.conjugate()*(Lambda12).transpose()*mH2I2) - 4*
      ms2*(Lambda12.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()*(
      Lambda12).transpose()) - 2*(Lambda12.conjugate()*(Lambda12).transpose()*
      TLambda12.conjugate()*(TLambda12).transpose()) - 2*(Lambda12.conjugate()*
      (TLambda12).transpose()*TLambda12.conjugate()*(Lambda12).transpose()) - 2
      *(TLambda12.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()*(
      TLambda12).transpose()) - 2*(TLambda12.conjugate()*(TLambda12).transpose(
      )*Lambda12.conjugate()*(Lambda12).transpose()) - mH2I2*Lambda12.conjugate
      ()*(Lambda12).transpose()*Lambda12.conjugate()*(Lambda12).transpose() - 2
      *(Lambda12.conjugate()*mH1I2.conjugate()*(Lambda12).transpose()*
      Lambda12.conjugate()*(Lambda12).transpose()) - 2*(Lambda12.conjugate()*(
      Lambda12).transpose()*mH2I2*Lambda12.conjugate()*(Lambda12).transpose())
      - 2*(Lambda12.conjugate()*(Lambda12).transpose()*Lambda12.conjugate()*
      mH1I2.conjugate()*(Lambda12).transpose()) - Lambda12.conjugate()*(
      Lambda12).transpose()*Lambda12.conjugate()*(Lambda12).transpose()*mH2I2 +
      0.04*(6*Conj(MassBp)*Sqr(gN)*(3*(MassB + 2*MassBp)*Sqr(g1) + 5*(2*MassBp
      + MassWB)*Sqr(g2) + 96*MassBp*Sqr(gN)) + 9*Conj(MassB)*Sqr(g1)*(99*MassB
      *Sqr(g1) + 5*(2*MassB + MassWB)*Sqr(g2) + 2*(2*MassB + MassBp)*Sqr(gN)) +
      5*(30*Power(g2,4)*Tr22 - 4.898979485566356*g1*gN*Tr2U114 -
      4.898979485566356*g1*gN*Tr2U141 + 15.491933384829668*g1*Tr31 -
      12.649110640673518*gN*Tr34 + 6*Tr2U111*Sqr(g1) + 4*Tr2U144*Sqr(gN) + 3*
      Conj(MassWB)*Sqr(g2)*(3*(MassB + 2*MassWB)*Sqr(g1) + 145*MassWB*Sqr(g2) +
      2*(MassBp + 2*MassWB)*Sqr(gN))))*UNITMATRIX(2))).real();


   return beta_mH2I2;
}

/**
 * Calculates the three-loop beta function of mH2I2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMtower_soft_parameters::calc_beta_mH2I2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,2,2> beta_mH2I2;

   beta_mH2I2 = ZEROMATRIX(2,2);


   return beta_mH2I2;
}

} // namespace flexiblesusy
