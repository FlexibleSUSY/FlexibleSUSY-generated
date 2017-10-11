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

// File generated at Tue 10 Oct 2017 21:25:29

#include "E6SSMEFTHiggs_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of mH1I2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMEFTHiggs_soft_parameters::calc_beta_mH1I2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,2,2> beta_mH1I2;

   beta_mH1I2 = (oneOver16PiSqr*(2*ms2*((Lambda12).adjoint()*Lambda12) +
      2*((TLambda12).adjoint()*TLambda12) + mH1I2*(Lambda12).adjoint()*Lambda12
      + 2*((Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) + (Lambda12)
      .adjoint()*Lambda12*mH1I2 - 0.1*(7.745966692414834*g1*Tr11 +
      9.486832980505138*gN*Tr14 + 12*AbsSqr(MassB)*Sqr(g1) + 60*AbsSqr(MassWB)*
      Sqr(g2) + 18*AbsSqr(MassBp)*Sqr(gN))*UNITMATRIX(2))).real();


   return beta_mH1I2;
}

/**
 * Calculates the 2-loop beta function of mH1I2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMEFTHiggs_soft_parameters::calc_beta_mH1I2_2_loop(const Soft_traces& soft_traces) const
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


   Eigen::Matrix<double,2,2> beta_mH1I2;

   beta_mH1I2 = (twoLoop*(2*(-3*traceconjTKappaTpTKappa - 2*
      traceconjTLambda12TpTLambda12 - 6*ms2*traceKappaAdjKappa - 3*
      traceKappaAdjKappaconjmDx2 - 3*traceKappaconjmDxbar2AdjKappa - 4*ms2*
      traceLambda12AdjLambda12 - 2*traceLambda12AdjLambda12conjmH2I2 - 2*
      tracemH1I2AdjLambda12Lambda12 - 2*(mHd2 + mHu2 + 2*ms2)*AbsSqr(Lambdax) -
      2*AbsSqr(TLambdax) + ms2*Sqr(gN) + 2*AbsSqr(MassBp)*Sqr(gN))*((Lambda12)
      .adjoint()*Lambda12) - 2*(3*traceconjTKappaTpKappa + 2*
      traceconjTLambda12TpLambda12 + 2*Conj(TLambdax)*Lambdax + Conj(MassBp)*
      Sqr(gN))*((Lambda12).adjoint()*TLambda12) - 2*(3*traceAdjKappaTKappa + 2*
      traceAdjLambda12TLambda12 + MassBp*Sqr(gN) + 2*Conj(Lambdax)*TLambdax)*((
      TLambda12).adjoint()*Lambda12) + 2*(-3*traceKappaAdjKappa - 2*
      traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + Sqr(gN))*((TLambda12)
      .adjoint()*TLambda12) + (-3*traceKappaAdjKappa - 2*
      traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + Sqr(gN))*(mH1I2*(Lambda12)
      .adjoint()*Lambda12) + 2*(-3*traceKappaAdjKappa - 2*
      traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + Sqr(gN))*((Lambda12)
      .adjoint()*mH2I2.conjugate()*Lambda12) + (-3*traceKappaAdjKappa - 2*
      traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + Sqr(gN))*((Lambda12)
      .adjoint()*Lambda12*mH1I2) - 4*ms2*((Lambda12).adjoint()*Lambda12*(
      Lambda12).adjoint()*Lambda12) - 2*((Lambda12).adjoint()*Lambda12*(
      TLambda12).adjoint()*TLambda12) - 2*((Lambda12).adjoint()*TLambda12*(
      TLambda12).adjoint()*Lambda12) - 2*((TLambda12).adjoint()*Lambda12*(
      Lambda12).adjoint()*TLambda12) - 2*((TLambda12).adjoint()*TLambda12*(
      Lambda12).adjoint()*Lambda12) - mH1I2*(Lambda12).adjoint()*Lambda12*(
      Lambda12).adjoint()*Lambda12 - 2*((Lambda12).adjoint()*mH2I2.conjugate()*
      Lambda12*(Lambda12).adjoint()*Lambda12) - 2*((Lambda12).adjoint()*
      Lambda12*mH1I2*(Lambda12).adjoint()*Lambda12) - 2*((Lambda12).adjoint()*
      Lambda12*(Lambda12).adjoint()*mH2I2.conjugate()*Lambda12) - (Lambda12)
      .adjoint()*Lambda12*(Lambda12).adjoint()*Lambda12*mH1I2 + 0.01*(9*Conj(
      MassBp)*Sqr(gN)*(-2*(MassB + 2*MassBp)*Sqr(g1) + 30*(2*MassBp + MassWB)*
      Sqr(g2) + 591*MassBp*Sqr(gN)) + 18*Conj(MassB)*Sqr(g1)*(198*MassB*Sqr(g1)
      + 10*(2*MassB + MassWB)*Sqr(g2) - (2*MassB + MassBp)*Sqr(gN)) + 10*(3*
      Conj(MassWB)*Sqr(g2)*(6*(MassB + 2*MassWB)*Sqr(g1) + 290*MassWB*Sqr(g2) +
      9*(MassBp + 2*MassWB)*Sqr(gN)) + 2*(7.348469228349534*g1*gN*Tr2U114 +
      7.348469228349534*g1*gN*Tr2U141 - 15.491933384829668*g1*Tr31 -
      18.973665961010276*gN*Tr34 + 30*Tr22*Quad(g2) + 6*Tr2U111*Sqr(g1) + 9*
      Tr2U144*Sqr(gN))))*UNITMATRIX(2))).real();


   return beta_mH1I2;
}

/**
 * Calculates the 3-loop beta function of mH1I2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,2,2> E6SSMEFTHiggs_soft_parameters::calc_beta_mH1I2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,2,2> beta_mH1I2;

   beta_mH1I2 = ZEROMATRIX(2,2);


   return beta_mH1I2;
}

} // namespace flexiblesusy