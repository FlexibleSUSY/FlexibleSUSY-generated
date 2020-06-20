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
 * Calculates the 1-loop beta function of me2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_me2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (4*mHd2*(Ye*Ye.adjoint()) + 4*(TYe*(TYe).adjoint()) + 2*(me2*Ye*
      Ye.adjoint()) + 4*(Ye*ml2*Ye.adjoint()) + 2*(Ye*Ye.adjoint()*me2) + 0.1*(
      15.491933384829668*g1*Tr11 + 3.1622776601683795*gN*Tr14 - 48*AbsSqr(MassB
      )*Sqr(g1) - 2*AbsSqr(MassBp)*Sqr(gN))*UNITMATRIX(3)).real();


   return oneLoop * beta_me2;
}

/**
 * Calculates the 2-loop beta function of me2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_me2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = (-0.8*(15*traceconjTYdTpTYd + 5*traceconjTYeTpTYe + 15*
      tracemd2YdAdjYd + 5*traceme2YeAdjYe + 5*traceml2AdjYeYe + 15*
      tracemq2AdjYdYd + 30*mHd2*traceYdAdjYd + 10*mHd2*traceYeAdjYe + 10*mHd2*
      AbsSqr(Lambdax) + 5*mHu2*AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*
      AbsSqr(TLambdax) + 3*mHd2*Sqr(g1) + 6*AbsSqr(MassB)*Sqr(g1) - 15*mHd2*Sqr
      (g2) - 30*AbsSqr(MassWB)*Sqr(g2) - 3*mHd2*Sqr(gN) - 6*AbsSqr(MassBp)*Sqr(
      gN))*(Ye*Ye.adjoint()) + 0.8*(-15*traceAdjYdTYd - 5*traceAdjYeTYe + 3*
      MassB*Sqr(g1) - 15*MassWB*Sqr(g2) - 3*MassBp*Sqr(gN) - 5*Conj(Lambdax)*
      TLambdax)*(Ye*(TYe).adjoint()) - 0.8*(15*traceconjTYdTpYd + 5*
      traceconjTYeTpYe + 5*Conj(TLambdax)*Lambdax - 3*Conj(MassB)*Sqr(g1) + 15*
      Conj(MassWB)*Sqr(g2) + 3*Conj(MassBp)*Sqr(gN))*(TYe*Ye.adjoint()) - 0.8*(
      15*traceYdAdjYd + 5*traceYeAdjYe + 5*AbsSqr(Lambdax) + 3*Sqr(g1) - 15*Sqr
      (g2) - 3*Sqr(gN))*(TYe*(TYe).adjoint()) - 0.4*(15*traceYdAdjYd + 5*
      traceYeAdjYe + 5*AbsSqr(Lambdax) + 3*Sqr(g1) - 15*Sqr(g2) - 3*Sqr(gN))*(
      me2*Ye*Ye.adjoint()) - 0.8*(15*traceYdAdjYd + 5*traceYeAdjYe + 5*AbsSqr(
      Lambdax) + 3*Sqr(g1) - 15*Sqr(g2) - 3*Sqr(gN))*(Ye*ml2*Ye.adjoint()) -
      0.4*(15*traceYdAdjYd + 5*traceYeAdjYe + 5*AbsSqr(Lambdax) + 3*Sqr(g1) -
      15*Sqr(g2) - 3*Sqr(gN))*(Ye*Ye.adjoint()*me2) - 8*mHd2*(Ye*Ye.adjoint()*
      Ye*Ye.adjoint()) - 4*(Ye*Ye.adjoint()*TYe*(TYe).adjoint()) - 4*(Ye*(TYe).
      adjoint()*TYe*Ye.adjoint()) - 4*(TYe*Ye.adjoint()*Ye*(TYe).adjoint()) - 4
      *(TYe*(TYe).adjoint()*Ye*Ye.adjoint()) - 2*(me2*Ye*Ye.adjoint()*Ye*Ye.
      adjoint()) - 4*(Ye*ml2*Ye.adjoint()*Ye*Ye.adjoint()) - 4*(Ye*Ye.adjoint()
      *me2*Ye*Ye.adjoint()) - 4*(Ye*Ye.adjoint()*Ye*ml2*Ye.adjoint()) - 2*(Ye*
      Ye.adjoint()*Ye*Ye.adjoint()*me2) + 0.01*(97.97958971132712*g1*gN*Tr2U114
       + 97.97958971132712*g1*gN*Tr2U141 + 619.6773353931867*g1*Tr31 +
      126.49110640673518*gN*Tr34 + 15552*AbsSqr(MassB)*Quad(g1) + 567*AbsSqr(
      MassBp)*Quad(gN) + 480*Tr2U111*Sqr(g1) + 20*Tr2U144*Sqr(gN) - 48*AbsSqr(
      MassB)*Sqr(g1)*Sqr(gN) - 48*AbsSqr(MassBp)*Sqr(g1)*Sqr(gN) - 24*MassBp*
      Conj(MassB)*Sqr(g1)*Sqr(gN) - 24*MassB*Conj(MassBp)*Sqr(g1)*Sqr(gN))*
      UNITMATRIX(3)).real();


   return twoLoop * beta_me2;
}

/**
 * Calculates the 3-loop beta function of me2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_me2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return threeLoop * beta_me2;
}

/**
 * Calculates the 4-loop beta function of me2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_me2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return fourLoop * beta_me2;
}

/**
 * Calculates the 5-loop beta function of me2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSM_soft_parameters::calc_beta_me2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_me2;

   beta_me2 = ZEROMATRIX(3,3);


   return fiveLoop * beta_me2;
}

} // namespace flexiblesusy
