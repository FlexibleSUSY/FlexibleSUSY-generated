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

// File generated at Mon 5 Mar 2018 16:04:56

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
 * Calculates the 1-loop beta function of mDxbar2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_mDxbar2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = (oneOver16PiSqr*(2*ms2*((Kappa).transpose()*
      Kappa.conjugate()) + 2*((TKappa).transpose()*TKappa.conjugate()) +
      mDxbar2*(Kappa).transpose()*Kappa.conjugate() + 2*((Kappa).transpose()*
      mDx2*Kappa.conjugate()) + (Kappa).transpose()*Kappa.conjugate()*mDxbar2 +
      0.03333333333333333*(15.491933384829668*g1*Tr11 - 28.460498941515414*gN*
      Tr14 - 16*AbsSqr(MassB)*Sqr(g1) - 320*AbsSqr(MassG)*Sqr(g3) - 54*AbsSqr(
      MassBp)*Sqr(gN))*UNITMATRIX(3))).real();


   return beta_mDxbar2;
}

/**
 * Calculates the 2-loop beta function of mDxbar2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_mDxbar2_2_loop(const Soft_traces& soft_traces) const
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
   const double traceAdjKappaTKappa = TRACE_STRUCT.traceAdjKappaTKappa;
   const double traceAdjLambda12TLambda12 =
      TRACE_STRUCT.traceAdjLambda12TLambda12;
   const double traceconjTKappaTpKappa =
      TRACE_STRUCT.traceconjTKappaTpKappa;
   const double traceconjTLambda12TpLambda12 =
      TRACE_STRUCT.traceconjTLambda12TpLambda12;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = (twoLoop*(2*(-3*traceconjTKappaTpTKappa - 2*
      traceconjTLambda12TpTLambda12 - 6*ms2*traceKappaAdjKappa - 3*
      traceKappaAdjKappaconjmDx2 - 3*traceKappaconjmDxbar2AdjKappa - 4*ms2*
      traceLambda12AdjLambda12 - 2*traceLambda12AdjLambda12conjmH2I2 - 2*
      tracemH1I2AdjLambda12Lambda12 - 2*(mHd2 + mHu2 + 2*ms2)*AbsSqr(Lambdax) -
      2*AbsSqr(TLambdax) + ms2*Sqr(gN) + 2*AbsSqr(MassBp)*Sqr(gN))*((Kappa)
      .transpose()*Kappa.conjugate()) - 2*(3*traceAdjKappaTKappa + 2*
      traceAdjLambda12TLambda12 + MassBp*Sqr(gN) + 2*Conj(Lambdax)*TLambdax)*((
      Kappa).transpose()*TKappa.conjugate()) - 2*(3*traceconjTKappaTpKappa + 2*
      traceconjTLambda12TpLambda12 + 2*Conj(TLambdax)*Lambdax + Conj(MassBp)*
      Sqr(gN))*((TKappa).transpose()*Kappa.conjugate()) + 2*(-3*
      traceKappaAdjKappa - 2*traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + Sqr
      (gN))*((TKappa).transpose()*TKappa.conjugate()) + (-3*traceKappaAdjKappa
      - 2*traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + Sqr(gN))*(mDxbar2*(
      Kappa).transpose()*Kappa.conjugate()) + 2*(-3*traceKappaAdjKappa - 2*
      traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + Sqr(gN))*((Kappa)
      .transpose()*mDx2*Kappa.conjugate()) + (-3*traceKappaAdjKappa - 2*
      traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + Sqr(gN))*((Kappa)
      .transpose()*Kappa.conjugate()*mDxbar2) - 4*ms2*((Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()) - 2*((Kappa)
      .transpose()*Kappa.conjugate()*(TKappa).transpose()*TKappa.conjugate()) -
      2*((Kappa).transpose()*TKappa.conjugate()*(TKappa).transpose()*
      Kappa.conjugate()) - 2*((TKappa).transpose()*Kappa.conjugate()*(Kappa)
      .transpose()*TKappa.conjugate()) - 2*((TKappa).transpose()*
      TKappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()) - mDxbar2*(
      Kappa).transpose()*Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate(
      ) - 2*((Kappa).transpose()*mDx2*Kappa.conjugate()*(Kappa).transpose()*
      Kappa.conjugate()) - 2*((Kappa).transpose()*Kappa.conjugate()*mDxbar2*(
      Kappa).transpose()*Kappa.conjugate()) - 2*((Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose()*mDx2*Kappa.conjugate()) - (Kappa)
      .transpose()*Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()*
      mDxbar2 + 0.0011111111111111111*(27*Conj(MassBp)*Sqr(gN)*(24*(MassB + 2*
      MassBp)*Sqr(g1) + 160*(2*MassBp + MassG)*Sqr(g3) + 1773*MassBp*Sqr(gN)) +
      8*Conj(MassB)*Sqr(g1)*(1752*MassB*Sqr(g1) + 160*(2*MassB + MassG)*Sqr(g3
      ) + 81*(2*MassB + MassBp)*Sqr(gN)) + 20*(3*(g1*(-14.696938456699067*gN*(
      Tr2U114 + Tr2U141) + 30.983866769659336*Tr31) + 9*gN*(3*gN*Tr2U144 -
      6.324555320336759*Tr34) + 160*Tr23*Quad(g3) + 8*Tr2U111*Sqr(g1)) + 8*Conj
      (MassG)*Sqr(g3)*(8*(MassB + 2*MassG)*Sqr(g1) + 300*MassG*Sqr(g3) + 27*(
      MassBp + 2*MassG)*Sqr(gN))))*UNITMATRIX(3))).real();


   return beta_mDxbar2;
}

/**
 * Calculates the 3-loop beta function of mDxbar2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_mDxbar2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = ZEROMATRIX(3,3);


   return beta_mDxbar2;
}

/**
 * Calculates the 4-loop beta function of mDxbar2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> CE6SSM_soft_parameters::calc_beta_mDxbar2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mDxbar2;

   beta_mDxbar2 = ZEROMATRIX(3,3);


   return beta_mDxbar2;
}

} // namespace flexiblesusy
