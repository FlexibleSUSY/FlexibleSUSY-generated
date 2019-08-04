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

// File generated at Sun 4 Aug 2019 17:19:08

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
 * Calculates the 1-loop beta function of md2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_md2_1_loop(const Soft_traces& soft_traces) const
{
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (oneOver16PiSqr*(4*mHd2*(Yd*Yd.adjoint()) + 4*(TYd*(TYd).adjoint(
      )) + 2*(md2*Yd*Yd.adjoint()) + 4*(Yd*mq2*Yd.adjoint()) + 2*(Yd*Yd.adjoint
      ()*md2) + 0.06666666666666667*(7.745966692414834*g1*Tr11 +
      9.486832980505138*gN*Tr14 - 8*AbsSqr(MassB)*Sqr(g1) - 160*AbsSqr(MassG)*
      Sqr(g3) - 12*AbsSqr(MassBp)*Sqr(gN))*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the 2-loop beta function of md2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_md2_2_loop(const Soft_traces& soft_traces) const
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

   beta_md2 = (twoLoop*(0.4*(-30*traceconjTYdTpTYd - 10*traceconjTYeTpTYe - 30*
      tracemd2YdAdjYd - 10*traceme2YeAdjYe - 10*traceml2AdjYeYe - 30*
      tracemq2AdjYdYd - 60*mHd2*traceYdAdjYd - 20*mHd2*traceYeAdjYe - 20*mHd2*
      AbsSqr(Lambdax) - 10*mHu2*AbsSqr(Lambdax) - 10*ms2*AbsSqr(Lambdax) - 10*
      AbsSqr(TLambdax) + 2*mHd2*Sqr(g1) + 4*AbsSqr(MassB)*Sqr(g1) + 30*mHd2*Sqr
      (g2) + 60*AbsSqr(MassWB)*Sqr(g2) + 3*mHd2*Sqr(gN) + 6*AbsSqr(MassBp)*Sqr(
      gN))*(Yd*Yd.adjoint()) - 0.4*(30*traceAdjYdTYd + 10*traceAdjYeTYe + 2*
      MassB*Sqr(g1) + 30*MassWB*Sqr(g2) + 3*MassBp*Sqr(gN) + 10*Conj(Lambdax)*
      TLambdax)*(Yd*(TYd).adjoint()) - 0.4*(30*traceconjTYdTpYd + 10*
      traceconjTYeTpYe + 10*Conj(TLambdax)*Lambdax + 2*Conj(MassB)*Sqr(g1) + 30
      *Conj(MassWB)*Sqr(g2) + 3*Conj(MassBp)*Sqr(gN))*(TYd*Yd.adjoint()) + 0.4*
      (-30*traceYdAdjYd - 10*traceYeAdjYe - 10*AbsSqr(Lambdax) + 2*Sqr(g1) + 30
      *Sqr(g2) + 3*Sqr(gN))*(TYd*(TYd).adjoint()) + 0.2*(-30*traceYdAdjYd - 10*
      traceYeAdjYe - 10*AbsSqr(Lambdax) + 2*Sqr(g1) + 30*Sqr(g2) + 3*Sqr(gN))*(
      md2*Yd*Yd.adjoint()) + 0.4*(-30*traceYdAdjYd - 10*traceYeAdjYe - 10*
      AbsSqr(Lambdax) + 2*Sqr(g1) + 30*Sqr(g2) + 3*Sqr(gN))*(Yd*mq2*Yd.adjoint(
      )) + 0.2*(-30*traceYdAdjYd - 10*traceYeAdjYe - 10*AbsSqr(Lambdax) + 2*Sqr
      (g1) + 30*Sqr(g2) + 3*Sqr(gN))*(Yd*Yd.adjoint()*md2) - 8*mHd2*(Yd*Yd.
      adjoint()*Yd*Yd.adjoint()) - 4*(Yd*Yd.adjoint()*TYd*(TYd).adjoint()) - 4*
      (mHd2 + mHu2)*(Yd*Yu.adjoint()*Yu*Yd.adjoint()) - 4*(Yd*Yu.adjoint()*TYu*
      (TYd).adjoint()) - 4*(Yd*(TYd).adjoint()*TYd*Yd.adjoint()) - 4*(Yd*(TYu).
      adjoint()*TYu*Yd.adjoint()) - 4*(TYd*Yd.adjoint()*Yd*(TYd).adjoint()) - 4
      *(TYd*Yu.adjoint()*Yu*(TYd).adjoint()) - 4*(TYd*(TYd).adjoint()*Yd*Yd.
      adjoint()) - 4*(TYd*(TYu).adjoint()*Yu*Yd.adjoint()) - 2*(md2*Yd*Yd.
      adjoint()*Yd*Yd.adjoint()) - 2*(md2*Yd*Yu.adjoint()*Yu*Yd.adjoint()) - 4*
      (Yd*mq2*Yd.adjoint()*Yd*Yd.adjoint()) - 4*(Yd*mq2*Yu.adjoint()*Yu*Yd.
      adjoint()) - 4*(Yd*Yd.adjoint()*md2*Yd*Yd.adjoint()) - 4*(Yd*Yd.adjoint()
      *Yd*mq2*Yd.adjoint()) - 2*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*md2) - 4*(Yd*
      Yu.adjoint()*mu2*Yu*Yd.adjoint()) - 4*(Yd*Yu.adjoint()*Yu*mq2*Yd.adjoint(
      )) - 2*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*md2) + 0.017777777777777778*(
      36.74234614174767*g1*gN*Tr2U114 + 36.74234614174767*g1*gN*Tr2U141 +
      116.1895003862225*g1*Tr31 + 142.30249470757707*gN*Tr34 + 876*AbsSqr(MassB
      )*Quad(g1) + 600*Tr23*Quad(g3) + 3000*AbsSqr(MassG)*Quad(g3) + 1296*
      AbsSqr(MassBp)*Quad(gN) + 30*Tr2U111*Sqr(g1) + 160*AbsSqr(MassB)*Sqr(g1)*
      Sqr(g3) + 160*AbsSqr(MassG)*Sqr(g1)*Sqr(g3) + 80*MassG*Conj(MassB)*Sqr(g1
      )*Sqr(g3) + 80*MassB*Conj(MassG)*Sqr(g1)*Sqr(g3) + 45*Tr2U144*Sqr(gN) -
      24*AbsSqr(MassB)*Sqr(g1)*Sqr(gN) - 24*AbsSqr(MassBp)*Sqr(g1)*Sqr(gN) - 12
      *MassBp*Conj(MassB)*Sqr(g1)*Sqr(gN) - 12*MassB*Conj(MassBp)*Sqr(g1)*Sqr(
      gN) + 240*AbsSqr(MassBp)*Sqr(g3)*Sqr(gN) + 240*AbsSqr(MassG)*Sqr(g3)*Sqr(
      gN) + 120*MassG*Conj(MassBp)*Sqr(g3)*Sqr(gN) + 120*MassBp*Conj(MassG)*Sqr
      (g3)*Sqr(gN))*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the 3-loop beta function of md2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_md2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

/**
 * Calculates the 4-loop beta function of md2.
 *
 * @return 4-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_md2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

/**
 * Calculates the 5-loop beta function of md2.
 *
 * @return 5-loop beta function
 */
Eigen::Matrix<double,3,3> E6SSMEFTHiggs_soft_parameters::calc_beta_md2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

} // namespace flexiblesusy
