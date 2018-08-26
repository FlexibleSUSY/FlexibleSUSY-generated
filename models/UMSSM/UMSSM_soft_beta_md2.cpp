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

// File generated at Sun 26 Aug 2018 14:23:15

#include "UMSSM_soft_parameters.hpp"
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
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = (oneOver16PiSqr*(4*mHd2*(Yd*Yd.adjoint()) + 4*(TYd*(TYd).adjoint(
      )) + 2*(md2*Yd*Yd.adjoint()) + 4*(Yd*mq2*Yd.adjoint()) + 2*(Yd*Yd.adjoint
      ()*md2) + 0.13333333333333333*(3.872983346207417*g1*Tr11 + 15*gp*Qd*Tr14
      - 4*AbsSqr(MassB)*Sqr(g1) - 80*AbsSqr(MassG)*Sqr(g3) - 60*AbsSqr(MassU)*
      Sqr(gp)*Sqr(Qd))*UNITMATRIX(3))).real();


   return beta_md2;
}

/**
 * Calculates the 2-loop beta function of md2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qs = INPUT(Qs);
   const auto Qu = INPUT(Qu);
   const auto Qv = INPUT(Qv);
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

   const Eigen::Matrix<double,3,3> beta_md2_1 = (UNITMATRIX(3)*(0.8*twoLoop*(-
      15*traceconjTYdTpTYd - 5*traceconjTYeTpTYe - 15*tracemd2YdAdjYd - 5*
      traceme2YeAdjYe - 5*traceml2AdjYeYe - 15*tracemq2AdjYdYd - 30*mHd2*
      traceYdAdjYd - 10*mHd2*traceYeAdjYe - 10*mHd2*AbsSqr(Lambdax) - 5*mHu2*
      AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) - 5*AbsSqr(TLambdax) + mHd2*Sqr(
      g1) + 2*AbsSqr(MassB)*Sqr(g1) + 15*mHd2*Sqr(g2) + 30*AbsSqr(MassWB)*Sqr(
      g2) - 10*mHd2*Sqr(gp)*Sqr(Qd) + 10*mHd2*Sqr(gp)*Sqr(QHd) - 20*AbsSqr(
      MassU)*Sqr(gp)*(Sqr(Qd) - Sqr(QHd) - Sqr(Qq)) + 10*mHd2*Sqr(gp)*Sqr(Qq))*
      (Yd*Yd.adjoint()) - 0.8*twoLoop*(MassB*Sqr(g1) + 5*(3*traceAdjYdTYd +
      traceAdjYeTYe + 3*MassWB*Sqr(g2) + 2*MassU*Sqr(gp)*(-Sqr(Qd) + Sqr(QHd) +
      Sqr(Qq))) + 5*Conj(Lambdax)*TLambdax)*(Yd*(TYd).adjoint()) - 0.8*twoLoop*
      (Conj(MassB)*Sqr(g1) + 5*(3*traceconjTYdTpYd + traceconjTYeTpYe + Conj(
      TLambdax)*Lambdax + 3*Conj(MassWB)*Sqr(g2) + 2*Conj(MassU)*Sqr(gp)*(-Sqr(
      Qd) + Sqr(QHd) + Sqr(Qq))))*(TYd*Yd.adjoint()) + 0.8*twoLoop*(-5*AbsSqr(
      Lambdax) + Sqr(g1) + 5*(-3*traceYdAdjYd - traceYeAdjYe + 3*Sqr(g2) + 2*
      Sqr(gp)*(-Sqr(Qd) + Sqr(QHd) + Sqr(Qq))))*(TYd*(TYd).adjoint()) + 0.4*
      twoLoop*(-5*AbsSqr(Lambdax) + Sqr(g1) + 5*(-3*traceYdAdjYd - traceYeAdjYe
       + 3*Sqr(g2) + 2*Sqr(gp)*(-Sqr(Qd) + Sqr(QHd) + Sqr(Qq))))*(md2*Yd*Yd.
      adjoint()) + 0.8*twoLoop*(-5*AbsSqr(Lambdax) + Sqr(g1) + 5*(-3*
      traceYdAdjYd - traceYeAdjYe + 3*Sqr(g2) + 2*Sqr(gp)*(-Sqr(Qd) + Sqr(QHd)
      + Sqr(Qq))))*(Yd*mq2*Yd.adjoint()) + 0.4*twoLoop*(-5*AbsSqr(Lambdax) +
      Sqr(g1) + 5*(-3*traceYdAdjYd - traceYeAdjYe + 3*Sqr(g2) + 2*Sqr(gp)*(-Sqr
      (Qd) + Sqr(QHd) + Sqr(Qq))))*(Yd*Yd.adjoint()*md2) - 8*mHd2*twoLoop*(Yd*
      Yd.adjoint()*Yd*Yd.adjoint()) - 4*twoLoop*(Yd*Yd.adjoint()*TYd*(TYd).
      adjoint()) - 4*(mHd2 + mHu2)*twoLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()) -
      4*twoLoop*(Yd*Yu.adjoint()*TYu*(TYd).adjoint()) - 4*twoLoop*(Yd*(TYd).
      adjoint()*TYd*Yd.adjoint()) - 4*twoLoop*(Yd*(TYu).adjoint()*TYu*Yd.
      adjoint()) - 4*twoLoop*(TYd*Yd.adjoint()*Yd*(TYd).adjoint()) - 4*twoLoop*
      (TYd*Yu.adjoint()*Yu*(TYd).adjoint()) - 4*twoLoop*(TYd*(TYd).adjoint()*Yd
      *Yd.adjoint()) - 4*twoLoop*(TYd*(TYu).adjoint()*Yu*Yd.adjoint()) - 2*
      twoLoop*(md2*Yd*Yd.adjoint()*Yd*Yd.adjoint()) - 2*twoLoop*(md2*Yd*Yu.
      adjoint()*Yu*Yd.adjoint()) - 4*twoLoop*(Yd*mq2*Yd.adjoint()*Yd*Yd.adjoint
      ()) - 4*twoLoop*(Yd*mq2*Yu.adjoint()*Yu*Yd.adjoint()) - 4*twoLoop*(Yd*Yd.
      adjoint()*md2*Yd*Yd.adjoint()) - 4*twoLoop*(Yd*Yd.adjoint()*Yd*mq2*Yd.
      adjoint()) - 2*twoLoop*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*md2) - 4*twoLoop*
      (Yd*Yu.adjoint()*mu2*Yu*Yd.adjoint()) - 4*twoLoop*(Yd*Yu.adjoint()*Yu*mq2
      *Yd.adjoint()) - 2*twoLoop*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*md2) +
      0.035555555555555556*twoLoop*(15*(3.872983346207417*g1*(gp*Qd*(Tr2U114 +
      Tr2U141) + Tr31) + 15*gp*Qd*(gp*Qd*Tr2U144 + Tr34) + 20*Tr23*Quad(g3) +
      Tr2U111*Sqr(g1)) + AbsSqr(MassB)*Sqr(g1)*(303*Sqr(g1) + 80*Sqr(g3)))*
      UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_md2_2 = (0.17777777777777778*twoLoop*(
      Conj(MassB)*Sqr(g1)*(8*MassG*Sqr(g3) + 3*(2*MassB + MassU)*Qd*(11*Qd + 3*
      (3*Qe - QHd + QHu - 3*Ql + 3*Qq - 6*Qu))*Sqr(gp)) + 8*Conj(MassG)*Sqr(g3)
      *((MassB + 2*MassG)*Sqr(g1) + 15*(-2*MassG*Sqr(g3) + (2*MassG + MassU)*
      Sqr(gp)*Sqr(Qd))) + 3*Qd*Conj(MassU)*Sqr(gp)*((MassB + 2*MassU)*(11*Qd +
      3*(3*Qe - QHd + QHu - 3*Ql + 3*Qq - 6*Qu))*Sqr(g1) + 5*Qd*(8*(MassG + 2*
      MassU)*Sqr(g3) + 9*MassU*Sqr(gp)*(11*Sqr(Qd) + 3*Sqr(Qe) + 2*Sqr(QHd) + 2
      *Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 9*Sqr(Qu) + 3*Sqr(Qv)))))*
      UNITMATRIX(3)).real();

   beta_md2 = beta_md2_1 + beta_md2_2;


   return beta_md2;
}

/**
 * Calculates the 3-loop beta function of md2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_3_loop(const Soft_traces& soft_traces) const
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
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_md2_4_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_md2;

   beta_md2 = ZEROMATRIX(3,3);


   return beta_md2;
}

} // namespace flexiblesusy
