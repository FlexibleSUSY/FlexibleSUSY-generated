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

// File generated at Sun 28 Aug 2016 15:10:26

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
 * Calculates the one-loop beta function of mDx2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_mDx2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mDx2;

   beta_mDx2 = (oneOver16PiSqr*(2*ms2*(Kappa.conjugate()*(Kappa)
      .transpose()) + 2*(TKappa.conjugate()*(TKappa).transpose()) + mDx2*
      Kappa.conjugate()*(Kappa).transpose() + 2*(Kappa.conjugate()*mDxbar2*(
      Kappa).transpose()) + Kappa.conjugate()*(Kappa).transpose()*mDx2 -
      0.5163977794943222*g1*Tr11*UNITMATRIX(3) - 0.6324555320336759*gN*Tr14*
      UNITMATRIX(3) - 0.5333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(3) -
      10.666666666666666*AbsSqr(MassG)*Sqr(g3)*UNITMATRIX(3) - 0.8*AbsSqr(
      MassBp)*Sqr(gN)*UNITMATRIX(3))).real();


   return beta_mDx2;
}

/**
 * Calculates the two-loop beta function of mDx2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_mDx2_two_loop(const Soft_traces& soft_traces) const
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
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mDx2;

   beta_mDx2 = (twoLoop*((-6*traceconjTKappaTpTKappa - 4*
      traceconjTLambda12TpTLambda12 - 12*ms2*traceKappaAdjKappa - 6*
      traceKappaAdjKappaconjmDx2 - 6*traceKappaconjmDxbar2AdjKappa - 8*ms2*
      traceLambda12AdjLambda12 - 4*traceLambda12AdjLambda12conjmH2I2 - 4*
      tracemH1I2AdjLambda12Lambda12 - 4*mHd2*AbsSqr(Lambdax) - 4*mHu2*AbsSqr(
      Lambdax) - 8*ms2*AbsSqr(Lambdax) - 4*AbsSqr(TLambdax) + 3*ms2*Sqr(gN) + 6
      *AbsSqr(MassBp)*Sqr(gN))*(Kappa.conjugate()*(Kappa).transpose()) + (-6*
      traceconjTKappaTpKappa - 4*traceconjTLambda12TpLambda12 - 4*Conj(TLambdax
      )*Lambdax - 3*Conj(MassBp)*Sqr(gN))*(Kappa.conjugate()*(TKappa).transpose
      ()) + (-6*traceAdjKappaTKappa - 4*traceAdjLambda12TLambda12 - 3*MassBp*
      Sqr(gN) - 4*Conj(Lambdax)*TLambdax)*(TKappa.conjugate()*(Kappa).transpose
      ()) + (-6*traceKappaAdjKappa - 4*traceLambda12AdjLambda12 - 4*AbsSqr(
      Lambdax) + 3*Sqr(gN))*(TKappa.conjugate()*(TKappa).transpose()) + (-3*
      traceKappaAdjKappa - 2*traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + 1.5
      *Sqr(gN))*(mDx2*Kappa.conjugate()*(Kappa).transpose()) + (-6*
      traceKappaAdjKappa - 4*traceLambda12AdjLambda12 - 4*AbsSqr(Lambdax) + 3*
      Sqr(gN))*(Kappa.conjugate()*mDxbar2*(Kappa).transpose()) + (-3*
      traceKappaAdjKappa - 2*traceLambda12AdjLambda12 - 2*AbsSqr(Lambdax) + 1.5
      *Sqr(gN))*(Kappa.conjugate()*(Kappa).transpose()*mDx2) - 4*ms2*(
      Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate()*(Kappa).transpose
      ()) - 2*(Kappa.conjugate()*(Kappa).transpose()*TKappa.conjugate()*(TKappa
      ).transpose()) - 2*(Kappa.conjugate()*(TKappa).transpose()*
      TKappa.conjugate()*(Kappa).transpose()) - 2*(TKappa.conjugate()*(Kappa)
      .transpose()*Kappa.conjugate()*(TKappa).transpose()) - 2*(
      TKappa.conjugate()*(TKappa).transpose()*Kappa.conjugate()*(Kappa)
      .transpose()) - mDx2*Kappa.conjugate()*(Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose() - 2*(Kappa.conjugate()*mDxbar2*(
      Kappa).transpose()*Kappa.conjugate()*(Kappa).transpose()) - 2*(
      Kappa.conjugate()*(Kappa).transpose()*mDx2*Kappa.conjugate()*(Kappa)
      .transpose()) - 2*(Kappa.conjugate()*(Kappa).transpose()*Kappa.conjugate(
      )*mDxbar2*(Kappa).transpose()) - Kappa.conjugate()*(Kappa).transpose()*
      Kappa.conjugate()*(Kappa).transpose()*mDx2 + 10.666666666666666*Power(g3,
      4)*Tr23*UNITMATRIX(3) + 0.6531972647421809*g1*gN*Tr2U114*UNITMATRIX(3) +
      0.6531972647421809*g1*gN*Tr2U141*UNITMATRIX(3) - 2.065591117977289*g1*
      Tr31*UNITMATRIX(3) - 2.5298221281347035*gN*Tr34*UNITMATRIX(3) +
      53.333333333333336*Power(g3,4)*AbsSqr(MassG)*UNITMATRIX(3) +
      0.5333333333333333*Tr2U111*Sqr(g1)*UNITMATRIX(3) + 2.8444444444444446*
      AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) + 1.4222222222222223*MassB*
      Conj(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) + 0.8*Tr2U144*Sqr(gN)*
      UNITMATRIX(3) + 4.266666666666667*AbsSqr(MassG)*Sqr(g3)*Sqr(gN)*
      UNITMATRIX(3) + 2.1333333333333333*MassBp*Conj(MassG)*Sqr(g3)*Sqr(gN)*
      UNITMATRIX(3) + 0.07111111111111111*Conj(MassB)*Sqr(g1)*(219*MassB*Sqr(g1
      ) + 20*(2*MassB + MassG)*Sqr(g3) - 3*(2*MassB + MassBp)*Sqr(gN))*
      UNITMATRIX(3) - 0.21333333333333335*Conj(MassBp)*Sqr(gN)*((MassB + 2*
      MassBp)*Sqr(g1) - 2*(5*(2*MassBp + MassG)*Sqr(g3) + 54*MassBp*Sqr(gN)))*
      UNITMATRIX(3))).real();


   return beta_mDx2;
}

/**
 * Calculates the three-loop beta function of mDx2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_mDx2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mDx2;

   beta_mDx2 = ZEROMATRIX(3,3);


   return beta_mDx2;
}

} // namespace flexiblesusy
