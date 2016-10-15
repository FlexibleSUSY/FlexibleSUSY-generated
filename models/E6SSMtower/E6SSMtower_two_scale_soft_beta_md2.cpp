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

// File generated at Sat 15 Oct 2016 16:01:28

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
 * Calculates the one-loop beta function of md2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMtower_soft_parameters::calc_beta_md2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (oneOver16PiSqr*(4*mHd2*(Yd*Yd.adjoint()) + 4*(TYd*(TYd)
      .adjoint()) + 2*(md2*Yd*Yd.adjoint()) + 4*(Yd*mq2*Yd.adjoint()) + 2*(Yd*
      Yd.adjoint()*md2) + 0.5163977794943222*g1*Tr11*UNITMATRIX(3) +
      0.6324555320336759*gN*Tr14*UNITMATRIX(3) - 0.5333333333333333*AbsSqr(
      MassB)*Sqr(g1)*UNITMATRIX(3) - 10.666666666666666*AbsSqr(MassG)*Sqr(g3)*
      UNITMATRIX(3) - 0.8*AbsSqr(MassBp)*Sqr(gN)*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the two-loop beta function of md2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMtower_soft_parameters::calc_beta_md2_two_loop(const Soft_traces& soft_traces) const
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
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (twoLoop*((-12*traceconjTYdTpTYd - 4*traceconjTYeTpTYe - 12
      *tracemd2YdAdjYd - 4*traceme2YeAdjYe - 4*traceml2AdjYeYe - 12*
      tracemq2AdjYdYd - 24*mHd2*traceYdAdjYd - 8*mHd2*traceYeAdjYe - 8*mHd2*
      AbsSqr(Lambdax) - 4*mHu2*AbsSqr(Lambdax) - 4*ms2*AbsSqr(Lambdax) - 4*
      AbsSqr(TLambdax) + 0.8*mHd2*Sqr(g1) + 1.6*AbsSqr(MassB)*Sqr(g1) + 12*mHd2
      *Sqr(g2) + 24*AbsSqr(MassWB)*Sqr(g2) + 1.2*mHd2*Sqr(gN) + 2.4*AbsSqr(
      MassBp)*Sqr(gN))*(Yd*Yd.adjoint()) + (-12*traceAdjYdTYd - 4*traceAdjYeTYe
      - 0.8*MassB*Sqr(g1) - 12*MassWB*Sqr(g2) - 1.2*MassBp*Sqr(gN) - 4*Conj(
      Lambdax)*TLambdax)*(Yd*(TYd).adjoint()) + (-12*traceconjTYdTpYd - 4*
      traceconjTYeTpYe - 4*Conj(TLambdax)*Lambdax - 0.8*Conj(MassB)*Sqr(g1) -
      12*Conj(MassWB)*Sqr(g2) - 1.2*Conj(MassBp)*Sqr(gN))*(TYd*Yd.adjoint()) +
      (-12*traceYdAdjYd - 4*traceYeAdjYe - 4*AbsSqr(Lambdax) + 0.8*Sqr(g1) + 12
      *Sqr(g2) + 1.2*Sqr(gN))*(TYd*(TYd).adjoint()) + (-6*traceYdAdjYd - 2*
      traceYeAdjYe - 2*AbsSqr(Lambdax) + 0.4*Sqr(g1) + 6*Sqr(g2) + 0.6*Sqr(gN))
      *(md2*Yd*Yd.adjoint()) + (-12*traceYdAdjYd - 4*traceYeAdjYe - 4*AbsSqr(
      Lambdax) + 0.8*Sqr(g1) + 12*Sqr(g2) + 1.2*Sqr(gN))*(Yd*mq2*Yd.adjoint())
      + (-6*traceYdAdjYd - 2*traceYeAdjYe - 2*AbsSqr(Lambdax) + 0.4*Sqr(g1) + 6
      *Sqr(g2) + 0.6*Sqr(gN))*(Yd*Yd.adjoint()*md2) - 8*mHd2*(Yd*Yd.adjoint()*
      Yd*Yd.adjoint()) - 4*(Yd*Yd.adjoint()*TYd*(TYd).adjoint()) + (-4*mHd2 - 4
      *mHu2)*(Yd*Yu.adjoint()*Yu*Yd.adjoint()) - 4*(Yd*Yu.adjoint()*TYu*(TYd)
      .adjoint()) - 4*(Yd*(TYd).adjoint()*TYd*Yd.adjoint()) - 4*(Yd*(TYu)
      .adjoint()*TYu*Yd.adjoint()) - 4*(TYd*Yd.adjoint()*Yd*(TYd).adjoint()) -
      4*(TYd*Yu.adjoint()*Yu*(TYd).adjoint()) - 4*(TYd*(TYd).adjoint()*Yd*
      Yd.adjoint()) - 4*(TYd*(TYu).adjoint()*Yu*Yd.adjoint()) - 2*(md2*Yd*
      Yd.adjoint()*Yd*Yd.adjoint()) - 2*(md2*Yd*Yu.adjoint()*Yu*Yd.adjoint()) -
      4*(Yd*mq2*Yd.adjoint()*Yd*Yd.adjoint()) - 4*(Yd*mq2*Yu.adjoint()*Yu*
      Yd.adjoint()) - 4*(Yd*Yd.adjoint()*md2*Yd*Yd.adjoint()) - 4*(Yd*
      Yd.adjoint()*Yd*mq2*Yd.adjoint()) - 2*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*
      md2) - 4*(Yd*Yu.adjoint()*mu2*Yu*Yd.adjoint()) - 4*(Yd*Yu.adjoint()*Yu*
      mq2*Yd.adjoint()) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*md2) +
      10.666666666666666*Power(g3,4)*Tr23*UNITMATRIX(3) + 0.6531972647421809*g1
      *gN*Tr2U114*UNITMATRIX(3) + 0.6531972647421809*g1*gN*Tr2U141*UNITMATRIX(3
      ) + 2.065591117977289*g1*Tr31*UNITMATRIX(3) + 2.5298221281347035*gN*Tr34*
      UNITMATRIX(3) + 53.333333333333336*Power(g3,4)*AbsSqr(MassG)*UNITMATRIX(3
      ) + 0.5333333333333333*Tr2U111*Sqr(g1)*UNITMATRIX(3) + 2.8444444444444446
      *AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) + 1.4222222222222223*MassB*
      Conj(MassG)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) + 0.8*Tr2U144*Sqr(gN)*
      UNITMATRIX(3) + 4.266666666666667*AbsSqr(MassG)*Sqr(g3)*Sqr(gN)*
      UNITMATRIX(3) + 2.1333333333333333*MassBp*Conj(MassG)*Sqr(g3)*Sqr(gN)*
      UNITMATRIX(3) + 0.07111111111111111*Conj(MassB)*Sqr(g1)*(219*MassB*Sqr(g1
      ) + 20*(2*MassB + MassG)*Sqr(g3) - 3*(2*MassB + MassBp)*Sqr(gN))*
      UNITMATRIX(3) - 0.21333333333333335*Conj(MassBp)*Sqr(gN)*((MassB + 2*
      MassBp)*Sqr(g1) - 2*(5*(2*MassBp + MassG)*Sqr(g3) + 54*MassBp*Sqr(gN)))*
      UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the three-loop beta function of md2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMtower_soft_parameters::calc_beta_md2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

} // namespace flexiblesusy
