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

// File generated at Fri 10 Apr 2020 20:18:23

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
 * Calculates the 1-loop beta function of mq2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mq2_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qq = INPUT(Qq);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = (oneOver16PiSqr*(2*mHd2*(Yd.adjoint()*Yd) + 2*mHu2*(Yu.adjoint()*
      Yu) + 2*((TYd).adjoint()*TYd) + 2*((TYu).adjoint()*TYu) + mq2*Yd.adjoint(
      )*Yd + mq2*Yu.adjoint()*Yu + 2*(Yd.adjoint()*md2*Yd) + Yd.adjoint()*Yd*
      mq2 + 2*(Yu.adjoint()*mu2*Yu) + Yu.adjoint()*Yu*mq2 + 0.06666666666666667
      *(3.872983346207417*g1*Tr11 + 30*gp*Qq*Tr14 - 2*AbsSqr(MassB)*Sqr(g1) -
      90*AbsSqr(MassWB)*Sqr(g2) - 160*AbsSqr(MassG)*Sqr(g3) - 120*AbsSqr(MassU)
      *Sqr(gp)*Sqr(Qq))*UNITMATRIX(3))).real();


   return beta_mq2;
}

/**
 * Calculates the 2-loop beta function of mq2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mq2_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qq = INPUT(Qq);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
   const auto Qv = INPUT(Qv);
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
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;
   const double traceAdjYdTYd = TRACE_STRUCT.traceAdjYdTYd;
   const double traceAdjYeTYe = TRACE_STRUCT.traceAdjYeTYe;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
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
      AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) + 2*mHd2*Sqr(g1) + 4*AbsSqr(MassB
      )*Sqr(g1) + 10*mHd2*Sqr(gp)*Sqr(Qd) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qd) +
      10*mHd2*Sqr(gp)*Sqr(QHd) + 20*AbsSqr(MassU)*Sqr(gp)*Sqr(QHd) - 10*mHd2*
      Sqr(gp)*Sqr(Qq) - 20*AbsSqr(MassU)*Sqr(gp)*Sqr(Qq))*(Yd.adjoint()*Yd) -
      0.4*twoLoop*(15*traceconjTYdTpYd + 5*traceconjTYeTpYe + 5*Conj(TLambdax)*
      Lambdax + 2*Conj(MassB)*Sqr(g1) + 10*Conj(MassU)*Sqr(gp)*Sqr(Qd) + 10*
      Conj(MassU)*Sqr(gp)*Sqr(QHd) - 10*Conj(MassU)*Sqr(gp)*Sqr(Qq))*(Yd.
      adjoint()*TYd) + 0.4*twoLoop*(-15*traceconjTYuTpTYu - 5*traceconjTYvTpTYv
       - 15*tracemq2AdjYuYu - 15*tracemu2YuAdjYu - 30*mHu2*traceYuAdjYu - 10*
      mHu2*traceYvAdjYv - 5*traceYvAdjYvconjml2 - 5*traceYvconjmvR2AdjYv - 5*
      mHd2*AbsSqr(Lambdax) - 10*mHu2*AbsSqr(Lambdax) - 5*ms2*AbsSqr(Lambdax) +
      4*mHu2*Sqr(g1) + 8*AbsSqr(MassB)*Sqr(g1) + 10*mHu2*Sqr(gp)*Sqr(QHu) + 20*
      AbsSqr(MassU)*Sqr(gp)*Sqr(QHu) - 10*mHu2*Sqr(gp)*Sqr(Qq) - 20*AbsSqr(
      MassU)*Sqr(gp)*Sqr(Qq) + 10*mHu2*Sqr(gp)*Sqr(Qu) + 20*AbsSqr(MassU)*Sqr(
      gp)*Sqr(Qu))*(Yu.adjoint()*Yu) - 0.4*twoLoop*(15*traceconjTYuTpYu + 5*
      traceconjTYvTpYv + 5*Conj(TLambdax)*Lambdax + 4*Conj(MassB)*Sqr(g1) + 10*
      Conj(MassU)*Sqr(gp)*Sqr(QHu) - 10*Conj(MassU)*Sqr(gp)*Sqr(Qq) + 10*Conj(
      MassU)*Sqr(gp)*Sqr(Qu))*(Yu.adjoint()*TYu) - 0.4*twoLoop*(15*
      traceAdjYdTYd + 5*traceAdjYeTYe + 2*MassB*Sqr(g1) + 10*MassU*Sqr(gp)*Sqr(
      Qd) + 10*MassU*Sqr(gp)*Sqr(QHd) - 10*MassU*Sqr(gp)*Sqr(Qq))*((TYd).
      adjoint()*Yd) + 0.4*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr
      (Lambdax) + 2*Sqr(g1) + 10*Sqr(gp)*Sqr(Qd) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr
      (gp)*Sqr(Qq))*((TYd).adjoint()*TYd) - 0.4*twoLoop*(15*traceAdjYuTYu + 5*
      traceAdjYvTYv + 4*MassB*Sqr(g1) + 10*MassU*Sqr(gp)*Sqr(QHu) - 10*MassU*
      Sqr(gp)*Sqr(Qq) + 10*MassU*Sqr(gp)*Sqr(Qu))*((TYu).adjoint()*Yu) + 0.4*
      twoLoop*(-15*traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) + 4*Sqr(g1
      ) + 10*Sqr(gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Qq) + 10*Sqr(gp)*Sqr(Qu))*((TYu)
      .adjoint()*TYu) + 0.2*twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*
      AbsSqr(Lambdax) + 2*Sqr(g1) + 10*Sqr(gp)*Sqr(Qd) + 10*Sqr(gp)*Sqr(QHd) -
      10*Sqr(gp)*Sqr(Qq))*(mq2*Yd.adjoint()*Yd) + 0.2*twoLoop*(-15*traceYuAdjYu
       - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) + 4*Sqr(g1) + 10*Sqr(gp)*Sqr(QHu) -
      10*Sqr(gp)*Sqr(Qq) + 10*Sqr(gp)*Sqr(Qu))*(mq2*Yu.adjoint()*Yu) + 0.4*
      twoLoop*(-15*traceYdAdjYd - 5*traceYeAdjYe - 5*AbsSqr(Lambdax) + 2*Sqr(g1
      ) + 10*Sqr(gp)*Sqr(Qd) + 10*Sqr(gp)*Sqr(QHd) - 10*Sqr(gp)*Sqr(Qq))*(Yd.
      adjoint()*md2*Yd) + 0.4*twoLoop*Sqr(g1)*(Yd.adjoint()*Yd*mq2))*UNITMATRIX
      (3)).real();
   const Eigen::Matrix<double,3,3> beta_mq2_2 = (UNITMATRIX(3)*(-2*twoLoop*
      AbsSqr(TLambdax)*(Yd.adjoint()*Yd) - 2*twoLoop*AbsSqr(TLambdax)*(Yu.
      adjoint()*Yu) - 2*twoLoop*Conj(Lambdax)*TLambdax*((TYd).adjoint()*Yd) - 2
      *twoLoop*Conj(Lambdax)*TLambdax*((TYu).adjoint()*Yu) + twoLoop*(-3*
      traceYdAdjYd - traceYeAdjYe - AbsSqr(Lambdax) + 2*Sqr(gp)*Sqr(Qd) + 2*Sqr
      (gp)*Sqr(QHd) - 2*Sqr(gp)*Sqr(Qq))*(Yd.adjoint()*Yd*mq2) + 0.4*twoLoop*(-
      15*traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(Lambdax) + 4*Sqr(g1) + 10*Sqr
      (gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Qq) + 10*Sqr(gp)*Sqr(Qu))*(Yu.adjoint()*
      mu2*Yu) + 0.2*twoLoop*(-15*traceYuAdjYu - 5*traceYvAdjYv - 5*AbsSqr(
      Lambdax) + 4*Sqr(g1) + 10*Sqr(gp)*Sqr(QHu) - 10*Sqr(gp)*Sqr(Qq) + 10*Sqr(
      gp)*Sqr(Qu))*(Yu.adjoint()*Yu*mq2) - 8*mHd2*twoLoop*(Yd.adjoint()*Yd*Yd.
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
      twoLoop*(Yu.adjoint()*Yu*Yu.adjoint()*Yu*mq2) + 0.0044444444444444444*
      twoLoop*(232.379000772445*g1*gp*Qq*Tr2U114 + 232.379000772445*g1*gp*Qq*
      Tr2U141 + 232.379000772445*g1*Tr31 + 1800*gp*Qq*Tr34 + 597*AbsSqr(MassB)*
      Quad(g1) + 1350*Tr22*Quad(g2) + 2400*Tr23*Quad(g3) - 9600*AbsSqr(MassG)*
      Quad(g3) + 108000*AbsSqr(MassU)*Quad(gp)*Quad(Qq) + 30*Tr2U111*Sqr(g1) +
      90*AbsSqr(MassB)*Sqr(g1)*Sqr(g2) + 45*MassWB*Conj(MassB)*Sqr(g1)*Sqr(g2)
      + 160*AbsSqr(MassB)*Sqr(g1)*Sqr(g3) + 160*AbsSqr(MassG)*Sqr(g1)*Sqr(g3) +
      80*MassG*Conj(MassB)*Sqr(g1)*Sqr(g3) + 80*MassB*Conj(MassG)*Sqr(g1)*Sqr(
      g3) + 7200*AbsSqr(MassG)*Sqr(g2)*Sqr(g3) + 3600*MassWB*Conj(MassG)*Sqr(g2
      )*Sqr(g3) + 1080*Qd*Qq*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 1080*Qe*Qq*AbsSqr(
      MassB)*Sqr(g1)*Sqr(gp) - 360*QHd*Qq*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 360*
      QHu*Qq*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) - 1080*Ql*Qq*AbsSqr(MassB)*Sqr(g1)*
      Sqr(gp) - 2160*Qq*Qu*AbsSqr(MassB)*Sqr(g1)*Sqr(gp) + 1080*Qd*Qq*AbsSqr(
      MassU)*Sqr(g1)*Sqr(gp) + 1080*Qe*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 360*
      QHd*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) + 360*QHu*Qq*AbsSqr(MassU)*Sqr(g1)*
      Sqr(gp) - 1080*Ql*Qq*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) + 540*MassU*Qd*Qq*Conj
      (MassB)*Sqr(g1)*Sqr(gp) + 540*MassU*Qe*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp) -
      180*MassU*QHd*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp) + 180*MassU*QHu*Qq*Conj(
      MassB)*Sqr(g1)*Sqr(gp) - 540*MassU*Ql*Qq*Conj(MassB)*Sqr(g1)*Sqr(gp) -
      1080*MassU*Qq*Qu*Conj(MassB)*Sqr(g1)*Sqr(gp) + 540*MassB*Qd*Qq*Conj(MassU
      )*Sqr(g1)*Sqr(gp) + 540*MassB*Qe*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp) - 180*
      MassB*QHd*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp) + 180*MassB*QHu*Qq*Conj(MassU)*
      Sqr(g1)*Sqr(gp) - 540*MassB*Ql*Qq*Conj(MassU)*Sqr(g1)*Sqr(gp) + 1800*
      Tr2U144*Sqr(gp)*Sqr(Qq) + 1200*AbsSqr(MassB)*Sqr(g1)*Sqr(gp)*Sqr(Qq) +
      1200*AbsSqr(MassU)*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 600*MassU*Conj(MassB)*Sqr(g1
      )*Sqr(gp)*Sqr(Qq) + 600*MassB*Conj(MassU)*Sqr(g1)*Sqr(gp)*Sqr(Qq) + 5400*
      AbsSqr(MassU)*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 2700*MassWB*Conj(MassU)*Sqr(g2)*
      Sqr(gp)*Sqr(Qq) + 9600*AbsSqr(MassG)*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 9600*
      AbsSqr(MassU)*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 4800*MassU*Conj(MassG)*Sqr(g3)*
      Sqr(gp)*Sqr(Qq) + 4800*MassG*Conj(MassU)*Sqr(g3)*Sqr(gp)*Sqr(Qq) + 48600*
      AbsSqr(MassU)*Quad(gp)*Sqr(Qd)*Sqr(Qq) + 16200*AbsSqr(MassU)*Quad(gp)*Sqr
      (Qe)*Sqr(Qq) + 10800*AbsSqr(MassU)*Quad(gp)*Sqr(QHd)*Sqr(Qq) + 10800*
      AbsSqr(MassU)*Quad(gp)*Sqr(QHu)*Sqr(Qq) + 32400*AbsSqr(MassU)*Quad(gp)*
      Sqr(Ql)*Sqr(Qq))*UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_mq2_3 = (0.2*twoLoop*(165*AbsSqr(MassWB
      )*Quad(g2) + 2*AbsSqr(MassWB)*Sqr(g1)*Sqr(g2) + MassB*Conj(MassWB)*Sqr(g1
      )*Sqr(g2) + 160*AbsSqr(MassWB)*Sqr(g2)*Sqr(g3) + 80*MassG*Conj(MassWB)*
      Sqr(g2)*Sqr(g3) - 48*Qq*Qu*AbsSqr(MassU)*Sqr(g1)*Sqr(gp) - 24*MassB*Qq*Qu
      *Conj(MassU)*Sqr(g1)*Sqr(gp) + 120*AbsSqr(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(Qq)
      + 60*MassU*Conj(MassWB)*Sqr(g2)*Sqr(gp)*Sqr(Qq) + 120*AbsSqr(MassU)*Quad(
      gp)*Sqr(Qq)*Sqr(Qs) + 1080*AbsSqr(MassU)*Quad(gp)*Sqr(Qq)*Sqr(Qu) + 360*
      AbsSqr(MassU)*Quad(gp)*Sqr(Qq)*Sqr(Qv))*UNITMATRIX(3)).real();

   beta_mq2 = beta_mq2_1 + beta_mq2_2 + beta_mq2_3;


   return beta_mq2;
}

/**
 * Calculates the 3-loop beta function of mq2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mq2_3_loop(const Soft_traces& soft_traces) const
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
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mq2_4_loop(const Soft_traces& soft_traces) const
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
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mq2_5_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mq2;

   beta_mq2 = ZEROMATRIX(3,3);


   return beta_mq2;
}

} // namespace flexiblesusy
