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

// File generated at Fri 10 Apr 2020 19:54:12

#include "TMSSM_soft_parameters.hpp"
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
 * Calculates the 1-loop beta function of ml2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_ml2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (oneOver16PiSqr*(2*mHd2*(Ye.adjoint()*Ye) + 2*((TYe).adjoint()*
      TYe) + ml2*Ye.adjoint()*Ye + 2*(Ye.adjoint()*me2*Ye) + Ye.adjoint()*Ye*
      ml2 - 0.2*(3.872983346207417*g1*Tr11 + 6*AbsSqr(MassB)*Sqr(g1) + 30*
      AbsSqr(MassWB)*Sqr(g2))*UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the 2-loop beta function of ml2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_ml2_2_loop(const Soft_traces& soft_traces) const
{
   const double traceYdAdjYd = TRACE_STRUCT.traceYdAdjYd;
   const double traceYeAdjYe = TRACE_STRUCT.traceYeAdjYe;
   const double traceconjTYdTpTYd = TRACE_STRUCT.traceconjTYdTpTYd;
   const double traceconjTYeTpTYe = TRACE_STRUCT.traceconjTYeTpTYe;
   const double tracemd2YdAdjYd = TRACE_STRUCT.tracemd2YdAdjYd;
   const double traceme2YeAdjYe = TRACE_STRUCT.traceme2YeAdjYe;
   const double traceml2AdjYeYe = TRACE_STRUCT.traceml2AdjYeYe;
   const double tracemq2AdjYdYd = TRACE_STRUCT.tracemq2AdjYdYd;
   const double traceconjTYdTpYd = TRACE_STRUCT.traceconjTYdTpYd;
   const double traceconjTYeTpYe = TRACE_STRUCT.traceconjTYeTpYe;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;


   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = (twoLoop*(0.2*(-30*traceconjTYdTpTYd - 10*traceconjTYeTpTYe - 30*
      tracemd2YdAdjYd - 10*traceme2YeAdjYe - 10*traceml2AdjYeYe - 30*
      tracemq2AdjYdYd - 60*mHd2*traceYdAdjYd - 20*mHd2*traceYeAdjYe - 30*mHd2*
      AbsSqr(Lambdax) - 15*mHu2*AbsSqr(Lambdax) - 15*mT2*AbsSqr(Lambdax) - 15*
      AbsSqr(TLambdax) + 12*mHd2*Sqr(g1) + 24*AbsSqr(MassB)*Sqr(g1))*(Ye.
      adjoint()*Ye) + 0.2*(-30*traceconjTYdTpYd - 10*traceconjTYeTpYe - 15*Conj
      (TLambdax)*Lambdax - 12*Conj(MassB)*Sqr(g1))*(Ye.adjoint()*TYe) + 0.2*(-
      30*traceAdjYdTYd - 10*traceAdjYeTYe - 12*MassB*Sqr(g1) - 15*Conj(Lambdax)
      *TLambdax)*((TYe).adjoint()*Ye) + 0.2*(-30*traceYdAdjYd - 10*traceYeAdjYe
       - 15*AbsSqr(Lambdax) + 12*Sqr(g1))*((TYe).adjoint()*TYe) + 0.1*(-30*
      traceYdAdjYd - 10*traceYeAdjYe - 15*AbsSqr(Lambdax) + 12*Sqr(g1))*(ml2*Ye
      .adjoint()*Ye) + 0.2*(-30*traceYdAdjYd - 10*traceYeAdjYe - 15*AbsSqr(
      Lambdax) + 12*Sqr(g1))*(Ye.adjoint()*me2*Ye) + 0.1*(-30*traceYdAdjYd - 10
      *traceYeAdjYe - 15*AbsSqr(Lambdax) + 12*Sqr(g1))*(Ye.adjoint()*Ye*ml2) -
      8*mHd2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 4*(Ye.adjoint()*Ye*(TYe).
      adjoint()*TYe) - 4*(Ye.adjoint()*TYe*(TYe).adjoint()*Ye) - 4*((TYe).
      adjoint()*Ye*Ye.adjoint()*TYe) - 4*((TYe).adjoint()*TYe*Ye.adjoint()*Ye)
      - 2*(ml2*Ye.adjoint()*Ye*Ye.adjoint()*Ye) - 4*(Ye.adjoint()*me2*Ye*Ye.
      adjoint()*Ye) - 4*(Ye.adjoint()*Ye*ml2*Ye.adjoint()*Ye) - 4*(Ye.adjoint()
      *Ye*Ye.adjoint()*me2*Ye) - 2*(Ye.adjoint()*Ye*Ye.adjoint()*Ye*ml2) + 0.04
      *(-77.45966692414834*g1*Tr31 + 621*AbsSqr(MassB)*Quad(g1) + 150*Tr22*Quad
      (g2) + 1725*AbsSqr(MassWB)*Quad(g2) + 30*Tr2U111*Sqr(g1) + 90*AbsSqr(
      MassB)*Sqr(g1)*Sqr(g2) + 90*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 45*MassWB*
      Conj(MassB)*Sqr(g1)*Sqr(g2) + 45*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2))*
      UNITMATRIX(3))).real();


   return beta_ml2;
}

/**
 * Calculates the 3-loop beta function of ml2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_ml2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

/**
 * Calculates the 4-loop beta function of ml2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_ml2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

/**
 * Calculates the 5-loop beta function of ml2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> TMSSM_soft_parameters::calc_beta_ml2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_ml2;

   beta_ml2 = ZEROMATRIX(3,3);


   return beta_ml2;
}

} // namespace flexiblesusy
