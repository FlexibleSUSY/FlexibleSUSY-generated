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

// File generated at Sun 28 Aug 2016 15:16:16

#include "NUTNMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mq2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> NUTNMSSM_soft_parameters::calc_beta_mq2_one_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;


   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = (oneOver16PiSqr*(2*mHd2*(Yd.adjoint()*Yd) + 2*mHu2*(
      Yu.adjoint()*Yu) + 2*((TYd).adjoint()*TYd) + 2*((TYu).adjoint()*TYu) +
      mq2*Yd.adjoint()*Yd + mq2*Yu.adjoint()*Yu + 2*(Yd.adjoint()*md2*Yd) +
      Yd.adjoint()*Yd*mq2 + 2*(Yu.adjoint()*mu2*Yu) + Yu.adjoint()*Yu*mq2 +
      0.2581988897471611*g1*Tr11*UNITMATRIX(3) - 0.13333333333333333*AbsSqr(
      MassB)*Sqr(g1)*UNITMATRIX(3) - 6*AbsSqr(MassWB)*Sqr(g2)*UNITMATRIX(3) -
      10.666666666666666*AbsSqr(MassG)*Sqr(g3)*UNITMATRIX(3))).real();


   return beta_mq2;
}

/**
 * Calculates the two-loop beta function of mq2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> NUTNMSSM_soft_parameters::calc_beta_mq2_two_loop(const Soft_traces& soft_traces) const
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
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr23 = TRACE_STRUCT.Tr23;


   Eigen::Matrix<double,3,3> beta_mq2;

   const Eigen::Matrix<double,3,3> beta_mq2_1 = (UNITMATRIX(3)*(0.4*
      twoLoop*(-15*traceconjTYdTpTYd - 5*traceconjTYeTpTYe - 15*tracemd2YdAdjYd
      - 5*traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*tracemq2AdjYdYd - 30*mHd2*
      traceYdAdjYd - 10*mHd2*traceYeAdjYe - 5*(2*mHd2 + mHu2 + ms2)*AbsSqr(
      Lambdax) - 5*AbsSqr(TLambdax) + 2*mHd2*Sqr(g1) + 4*AbsSqr(MassB)*Sqr(g1))
      *(Yd.adjoint()*Yd) - 0.4*twoLoop*(5*(3*traceconjTYdTpYd +
      traceconjTYeTpYe + Conj(TLambdax)*Lambdax) + 2*Conj(MassB)*Sqr(g1))*(
      Yd.adjoint()*TYd) + twoLoop*(-6*traceconjTYuTpTYu - 6*tracemq2AdjYuYu - 6
      *tracemu2YuAdjYu - 12*mHu2*traceYuAdjYu - 2*mHd2*AbsSqr(Lambdax) - 4*mHu2
      *AbsSqr(Lambdax) - 2*ms2*AbsSqr(Lambdax) - 2*AbsSqr(TLambdax) + 1.6*mHu2*
      Sqr(g1) + 3.2*AbsSqr(MassB)*Sqr(g1))*(Yu.adjoint()*Yu) + twoLoop*(-6*
      traceconjTYuTpYu - 2*Conj(TLambdax)*Lambdax - 1.6*Conj(MassB)*Sqr(g1))*(
      Yu.adjoint()*TYu) + twoLoop*(-6*traceAdjYdTYd - 2*traceAdjYeTYe - 0.8*
      MassB*Sqr(g1) - 2*Conj(Lambdax)*TLambdax)*((TYd).adjoint()*Yd) + twoLoop*
      (-6*traceYdAdjYd - 2*traceYeAdjYe - 2*AbsSqr(Lambdax) + 0.8*Sqr(g1))*((
      TYd).adjoint()*TYd) + twoLoop*(-6*traceAdjYuTYu - 1.6*MassB*Sqr(g1) - 2*
      Conj(Lambdax)*TLambdax)*((TYu).adjoint()*Yu) + twoLoop*(-6*traceYuAdjYu -
      2*AbsSqr(Lambdax) + 1.6*Sqr(g1))*((TYu).adjoint()*TYu) + twoLoop*(-3*
      traceYdAdjYd - traceYeAdjYe - AbsSqr(Lambdax) + 0.4*Sqr(g1))*(mq2*
      Yd.adjoint()*Yd) + twoLoop*(-3*traceYuAdjYu - AbsSqr(Lambdax) + 0.8*Sqr(
      g1))*(mq2*Yu.adjoint()*Yu) + twoLoop*(-6*traceYdAdjYd - 2*traceYeAdjYe -
      2*AbsSqr(Lambdax) + 0.8*Sqr(g1))*(Yd.adjoint()*md2*Yd) + twoLoop*(-3*
      traceYdAdjYd - traceYeAdjYe - AbsSqr(Lambdax) + 0.4*Sqr(g1))*(Yd.adjoint(
      )*Yd*mq2) + twoLoop*(-6*traceYuAdjYu - 2*AbsSqr(Lambdax) + 1.6*Sqr(g1))*(
      Yu.adjoint()*mu2*Yu) + twoLoop*(-3*traceYuAdjYu - AbsSqr(Lambdax) + 0.8*
      Sqr(g1))*(Yu.adjoint()*Yu*mq2) - 8*mHd2*twoLoop*(Yd.adjoint()*Yd*
      Yd.adjoint()*Yd) - 4*twoLoop*(Yd.adjoint()*Yd*(TYd).adjoint()*TYd) - 4*
      twoLoop*(Yd.adjoint()*TYd*(TYd).adjoint()*Yd) - 8*mHu2*twoLoop*(
      Yu.adjoint()*Yu*Yu.adjoint()*Yu) - 4*twoLoop*(Yu.adjoint()*Yu*(TYu)
      .adjoint()*TYu) - 4*twoLoop*(Yu.adjoint()*TYu*(TYu).adjoint()*Yu) - 4*
      twoLoop*((TYd).adjoint()*Yd*Yd.adjoint()*TYd) - 4*twoLoop*((TYd).adjoint(
      )*TYd*Yd.adjoint()*Yd) - 4*twoLoop*((TYu).adjoint()*Yu*Yu.adjoint()*TYu)
      - 4*twoLoop*((TYu).adjoint()*TYu*Yu.adjoint()*Yu) - 2*twoLoop*(mq2*
      Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 2*twoLoop*(mq2*Yu.adjoint()*Yu*
      Yu.adjoint()*Yu) - 4*twoLoop*(Yd.adjoint()*md2*Yd*Yd.adjoint()*Yd) - 4*
      twoLoop*(Yd.adjoint()*Yd*mq2*Yd.adjoint()*Yd) - 4*twoLoop*(Yd.adjoint()*
      Yd*Yd.adjoint()*md2*Yd) - 2*twoLoop*(Yd.adjoint()*Yd*Yd.adjoint()*Yd*mq2)
      - 4*twoLoop*(Yu.adjoint()*mu2*Yu*Yu.adjoint()*Yu) - 4*twoLoop*(
      Yu.adjoint()*Yu*mq2*Yu.adjoint()*Yu) - 4*twoLoop*(Yu.adjoint()*Yu*
      Yu.adjoint()*mu2*Yu) - 2*twoLoop*(Yu.adjoint()*Yu*Yu.adjoint()*Yu*mq2) +
      twoLoop*(6*Power(g2,4)*Tr22*UNITMATRIX(3) + 10.666666666666666*Power(g3,4
      )*Tr23*UNITMATRIX(3) + 1.0327955589886444*g1*Tr31*UNITMATRIX(3) +
      2.6533333333333333*Power(g1,4)*AbsSqr(MassB)*UNITMATRIX(3) -
      42.666666666666664*Power(g3,4)*AbsSqr(MassG)*UNITMATRIX(3) +
      0.13333333333333333*Tr2U111*Sqr(g1)*UNITMATRIX(3) + 0.4*AbsSqr(MassB)*Sqr
      (g1)*Sqr(g2)*UNITMATRIX(3) + 0.2*MassWB*Conj(MassB)*Sqr(g1)*Sqr(g2)*
      UNITMATRIX(3) + 0.7111111111111111*AbsSqr(MassB)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 0.7111111111111111*AbsSqr(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 0.35555555555555557*MassG*Conj(MassB)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 0.35555555555555557*MassB*Conj(MassG)*Sqr(g1)*Sqr(g3)*
      UNITMATRIX(3) + 32*AbsSqr(MassG)*Sqr(g2)*Sqr(g3)*UNITMATRIX(3) + 16*
      MassWB*Conj(MassG)*Sqr(g2)*Sqr(g3)*UNITMATRIX(3)))).real();
   const Eigen::Matrix<double,3,3> beta_mq2_2 = (0.2*twoLoop*Conj(MassWB)
      *Sqr(g2)*((MassB + 2*MassWB)*Sqr(g1) + 5*(33*MassWB*Sqr(g2) + 16*(MassG +
      2*MassWB)*Sqr(g3)))*UNITMATRIX(3)).real();

   beta_mq2 = beta_mq2_1 + beta_mq2_2;


   return beta_mq2;
}

/**
 * Calculates the three-loop beta function of mq2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> NUTNMSSM_soft_parameters::calc_beta_mq2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = ZEROMATRIX(3,3);


   return beta_mq2;
}

} // namespace flexiblesusy
