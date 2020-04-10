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

// File generated at Fri 10 Apr 2020 18:30:49

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
 * Calculates the 1-loop beta function of mq2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_mq2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = (oneOver16PiSqr*(2*mHd2*(Yd.adjoint()*Yd) + 2*mHu2*(Yu.adjoint()*
      Yu) + 2*((TYd).adjoint()*TYd) + 2*((TYu).adjoint()*TYu) + mq2*Yd.adjoint(
      )*Yd + mq2*Yu.adjoint()*Yu + 2*(Yd.adjoint()*md2*Yd) + Yd.adjoint()*Yd*
      mq2 + 2*(Yu.adjoint()*mu2*Yu) + Yu.adjoint()*Yu*mq2 + 0.03333333333333333
      *(7.745966692414834*g1*Tr11 + 9.486832980505138*gN*Tr14 - 4*AbsSqr(MassB)
      *Sqr(g1) - 180*AbsSqr(MassWB)*Sqr(g2) - 320*AbsSqr(MassG)*Sqr(g3) - 6*
      AbsSqr(MassBp)*Sqr(gN))*UNITMATRIX(3))).real();


   return beta_mq2;
}

/**
 * Calculates the 2-loop beta function of mq2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_mq2_2_loop(const Soft_traces& soft_traces) const
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
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr22 = TRACE_STRUCT.Tr22;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mq2;

   const Eigen::Matrix<double,3,3> beta_mq2_1 = ((0.4*twoLoop*(-15*
      traceconjTYdTpTYd - 5*traceconjTYeTpTYe - 15*tracemd2YdAdjYd - 5*
      traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*tracemq2AdjYdYd - 30*mHd2*
      traceYdAdjYd - 10*mHd2*traceYeAdjYe - 10*mHd2*AbsSqr(Lambdax) - 5*mHu2*
      AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) - 5*AbsSqr(TLambdax) + 2*mHd2*Sqr
      (g1) + 4*AbsSqr(MassB)*Sqr(g1) + 3*mHd2*Sqr(gN) + 6*AbsSqr(MassBp)*Sqr(gN
      ))*(Yd.adjoint()*Yd) - 0.4*twoLoop*(15*traceconjTYdTpYd + 5*
      traceconjTYeTpYe + 5*Conj(TLambdax)*Lambdax + 2*Conj(MassB)*Sqr(g1) + 3*
      Conj(MassBp)*Sqr(gN))*(Yd.adjoint()*TYd) + 0.4*twoLoop*(-15*
      traceconjTYuTpTYu - 15*tracemq2AdjYuYu - 15*tracemu2YuAdjYu - 30*mHu2*
      traceYuAdjYu - 5*mHd2*AbsSqr(Lambdax) - 10*mHu2*AbsSqr(Lambdax) - 5*ms2*
      AbsSqr(Lambdax) - 5*AbsSqr(TLambdax) + 4*mHu2*Sqr(g1) + 8*AbsSqr(MassB)*
      Sqr(g1) + mHu2*Sqr(gN) + 2*AbsSqr(MassBp)*Sqr(gN))*(Yu.adjoint()*Yu) -
      0.4*twoLoop*(15*traceconjTYuTpYu + 5*Conj(TLambdax)*Lambdax + 4*Conj(
      MassB)*Sqr(g1) + Conj(MassBp)*Sqr(gN))*(Yu.adjoint()*TYu) - 0.4*twoLoop*(
      15*traceAdjYdTYd + 5*traceAdjYeTYe + 2*MassB*Sqr(g1) + 3*MassBp*Sqr(gN))*
      ((TYd).adjoint()*Yd) + 0.4*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5
      *AbsSqr(Lambdax) + 2*Sqr(g1) + 3*Sqr(gN))*((TYd).adjoint()*TYd) - 0.4*
      twoLoop*(15*traceAdjYuTYu + 4*MassB*Sqr(g1) + MassBp*Sqr(gN))*((TYu).
      adjoint()*Yu) + 0.4*twoLoop*(-15*traceYuAdjYu - 5*AbsSqr(Lambdax) + 4*Sqr
      (g1) + Sqr(gN))*((TYu).adjoint()*TYu) + 0.2*twoLoop*(-15*traceYdAdjYd - 5
      *traceYeAdjYe - 5*AbsSqr(Lambdax) + 2*Sqr(g1) + 3*Sqr(gN))*(mq2*Yd.
      adjoint()*Yd) + 0.2*twoLoop*(-15*traceYuAdjYu - 5*AbsSqr(Lambdax) + 4*Sqr
      (g1) + Sqr(gN))*(mq2*Yu.adjoint()*Yu) + 0.4*twoLoop*(-15*traceYdAdjYd - 5
      *traceYeAdjYe - 5*AbsSqr(Lambdax) + 2*Sqr(g1) + 3*Sqr(gN))*(Yd.adjoint()*
      md2*Yd) + 0.2*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(
      Lambdax) + 2*Sqr(g1) + 3*Sqr(gN))*(Yd.adjoint()*Yd*mq2) + 0.4*twoLoop*(-
      15*traceYuAdjYu - 5*AbsSqr(Lambdax) + 4*Sqr(g1) + Sqr(gN))*(Yu.adjoint()*
      mu2*Yu) + 0.2*twoLoop*(-15*traceYuAdjYu - 5*AbsSqr(Lambdax) + 4*Sqr(g1) +
      Sqr(gN))*(Yu.adjoint()*Yu*mq2) - 8*mHd2*twoLoop*(Yd.adjoint()*Yd*Yd.
      adjoint()*Yd) - 4*twoLoop*(Yd.adjoint()*Yd*(TYd).adjoint()*TYd) - 4*
      twoLoop*(Yd.adjoint()*TYd*(TYd).adjoint()*Yd) - 8*mHu2*twoLoop*(Yu.
      adjoint()*Yu*Yu.adjoint()*Yu) - 4*twoLoop*(Yu.adjoint()*Yu*(TYu).adjoint(
      )*TYu) - 4*twoLoop*(Yu.adjoint()*TYu*(TYu).adjoint()*Yu) - 4*twoLoop*((
      TYd).adjoint()*Yd*Yd.adjoint()*TYd) - 4*twoLoop*((TYd).adjoint()*TYd*Yd.
      adjoint()*Yd) - 4*twoLoop*((TYu).adjoint()*Yu*Yu.adjoint()*TYu) - 4*
      twoLoop*((TYu).adjoint()*TYu*Yu.adjoint()*Yu) - 2*twoLoop*(mq2*Yd.adjoint
      ()*Yd*Yd.adjoint()*Yd) - 2*twoLoop*(mq2*Yu.adjoint()*Yu*Yu.adjoint()*Yu)
      - 4*twoLoop*(Yd.adjoint()*md2*Yd*Yd.adjoint()*Yd) - 4*twoLoop*(Yd.adjoint
      ()*Yd*mq2*Yd.adjoint()*Yd) - 4*twoLoop*(Yd.adjoint()*Yd*Yd.adjoint()*md2*
      Yd) - 2*twoLoop*(Yd.adjoint()*Yd*Yd.adjoint()*Yd*mq2) - 4*twoLoop*(Yu.
      adjoint()*mu2*Yu*Yu.adjoint()*Yu) - 4*twoLoop*(Yu.adjoint()*Yu*mq2*Yu.
      adjoint()*Yu) - 4*twoLoop*(Yu.adjoint()*Yu*Yu.adjoint()*mu2*Yu) - 2*
      twoLoop*(Yu.adjoint()*Yu*Yu.adjoint()*Yu*mq2))*UNITMATRIX(3)).real();
   const Eigen::Matrix<double,3,3> beta_mq2_2 = (UNITMATRIX(3)*(-2*twoLoop*Conj
      (Lambdax)*TLambdax*((TYd).adjoint()*Yd) - 2*twoLoop*Conj(Lambdax)*
      TLambdax*((TYu).adjoint()*Yu) + 0.0011111111111111111*twoLoop*(
      146.96938456699067*g1*gN*Tr2U114 + 146.96938456699067*g1*gN*Tr2U141 +
      929.51600308978*g1*Tr31 + 1138.4199576606165*gN*Tr34 + 3468*AbsSqr(MassB)
      *Quad(g1) + 5400*Tr22*Quad(g2) + 78300*AbsSqr(MassWB)*Quad(g2) + 9600*
      Tr23*Quad(g3) + 48000*AbsSqr(MassG)*Quad(g3) + 5103*AbsSqr(MassBp)*Quad(
      gN) + 120*Tr2U111*Sqr(g1) + 360*AbsSqr(MassB)*Sqr(g1)*Sqr(g2) + 360*
      AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + 180*MassWB*Conj(MassB)*Sqr(g1)*Sqr(g2) +
      180*MassB*Conj(MassWB)*Sqr(g1)*Sqr(g2) + 640*AbsSqr(MassB)*Sqr(g1)*Sqr(g3
      ) + 640*AbsSqr(MassG)*Sqr(g1)*Sqr(g3) + 320*MassG*Conj(MassB)*Sqr(g1)*Sqr
      (g3) + 320*MassB*Conj(MassG)*Sqr(g1)*Sqr(g3) + 28800*AbsSqr(MassG)*Sqr(g2
      )*Sqr(g3) + 28800*AbsSqr(MassWB)*Sqr(g2)*Sqr(g3) + 14400*MassWB*Conj(
      MassG)*Sqr(g2)*Sqr(g3) + 14400*MassG*Conj(MassWB)*Sqr(g2)*Sqr(g3) + 180*
      Tr2U144*Sqr(gN) - 132*AbsSqr(MassB)*Sqr(g1)*Sqr(gN) - 132*AbsSqr(MassBp)*
      Sqr(g1)*Sqr(gN) - 66*MassBp*Conj(MassB)*Sqr(g1)*Sqr(gN) - 66*MassB*Conj(
      MassBp)*Sqr(g1)*Sqr(gN) + 540*AbsSqr(MassBp)*Sqr(g2)*Sqr(gN) + 540*AbsSqr
      (MassWB)*Sqr(g2)*Sqr(gN) + 270*MassWB*Conj(MassBp)*Sqr(g2)*Sqr(gN) + 270*
      MassBp*Conj(MassWB)*Sqr(g2)*Sqr(gN) + 960*AbsSqr(MassBp)*Sqr(g3)*Sqr(gN)
      + 960*AbsSqr(MassG)*Sqr(g3)*Sqr(gN) + 480*MassG*Conj(MassBp)*Sqr(g3)*Sqr(
      gN) + 480*MassBp*Conj(MassG)*Sqr(g3)*Sqr(gN))*UNITMATRIX(3))).real();

   beta_mq2 = beta_mq2_1 + beta_mq2_2;


   return beta_mq2;
}

/**
 * Calculates the 3-loop beta function of mq2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_mq2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = ZEROMATRIX(3,3);


   return beta_mq2;
}

/**
 * Calculates the 4-loop beta function of mq2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_mq2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = ZEROMATRIX(3,3);


   return beta_mq2;
}

/**
 * Calculates the 5-loop beta function of mq2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_mq2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = ZEROMATRIX(3,3);


   return beta_mq2;
}

} // namespace flexiblesusy
