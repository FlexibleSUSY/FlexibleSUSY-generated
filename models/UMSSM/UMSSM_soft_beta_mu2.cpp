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

// File generated at Tue 10 Oct 2017 22:15:55

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
 * Calculates the 1-loop beta function of mu2.
 *
 * @return 1-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_1_loop(const Soft_traces& soft_traces) const
{
   const auto Qu = INPUT(Qu);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (oneOver16PiSqr*(4*mHu2*(Yu*Yu.adjoint()) + 4*(TYu*(TYu)
      .adjoint()) + 2*(mu2*Yu*Yu.adjoint()) + 4*(Yu*mq2*Yu.adjoint()) + 2*(Yu*
      Yu.adjoint()*mu2) - 0.13333333333333333*(7.745966692414834*g1*Tr11 - 15*
      gp*Qu*Tr14 + 16*AbsSqr(MassB)*Sqr(g1) + 80*AbsSqr(MassG)*Sqr(g3) + 60*
      AbsSqr(MassU)*Sqr(gp)*Sqr(Qu))*UNITMATRIX(3))).real();


   return beta_mu2;
}

/**
 * Calculates the 2-loop beta function of mu2.
 *
 * @return 2-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_2_loop(const Soft_traces& soft_traces) const
{
   const auto Qd = INPUT(Qd);
   const auto Qe = INPUT(Qe);
   const auto QHd = INPUT(QHd);
   const auto QHu = INPUT(QHu);
   const auto Ql = INPUT(Ql);
   const auto Qq = INPUT(Qq);
   const auto Qu = INPUT(Qu);
   const auto Qs = INPUT(Qs);
   const auto Qv = INPUT(Qv);
   const double traceYuAdjYu = TRACE_STRUCT.traceYuAdjYu;
   const double traceYvAdjYv = TRACE_STRUCT.traceYvAdjYv;
   const double traceconjTYuTpTYu = TRACE_STRUCT.traceconjTYuTpTYu;
   const double traceconjTYvTpTYv = TRACE_STRUCT.traceconjTYvTpTYv;
   const double tracemq2AdjYuYu = TRACE_STRUCT.tracemq2AdjYuYu;
   const double tracemu2YuAdjYu = TRACE_STRUCT.tracemu2YuAdjYu;
   const double traceYvAdjYvconjml2 = TRACE_STRUCT.traceYvAdjYvconjml2;
   const double traceYvconjmvR2AdjYv = TRACE_STRUCT.traceYvconjmvR2AdjYv;
   const double traceAdjYuTYu = TRACE_STRUCT.traceAdjYuTYu;
   const double traceAdjYvTYv = TRACE_STRUCT.traceAdjYvTYv;
   const double traceconjTYuTpYu = TRACE_STRUCT.traceconjTYuTpYu;
   const double traceconjTYvTpYv = TRACE_STRUCT.traceconjTYvTpYv;
   const double Tr2U111 = TRACE_STRUCT.Tr2U111;
   const double Tr2U114 = TRACE_STRUCT.Tr2U114;
   const double Tr31 = TRACE_STRUCT.Tr31;
   const double Tr23 = TRACE_STRUCT.Tr23;
   const double Tr2U141 = TRACE_STRUCT.Tr2U141;
   const double Tr2U144 = TRACE_STRUCT.Tr2U144;
   const double Tr34 = TRACE_STRUCT.Tr34;


   Eigen::Matrix<double,3,3> beta_mu2;

   const Eigen::Matrix<double,3,3> beta_mu2_1 = (UNITMATRIX(3)*(-0.8*
      twoLoop*(15*traceconjTYuTpTYu + 5*traceconjTYvTpTYv + 15*tracemq2AdjYuYu
      + 15*tracemu2YuAdjYu + 30*mHu2*traceYuAdjYu + 10*mHu2*traceYvAdjYv + 5*
      traceYvAdjYvconjml2 + 5*traceYvconjmvR2AdjYv + 5*mHd2*AbsSqr(Lambdax) +
      10*mHu2*AbsSqr(Lambdax) + 5*ms2*AbsSqr(Lambdax) + 5*AbsSqr(TLambdax) +
      mHu2*Sqr(g1) + 2*AbsSqr(MassB)*Sqr(g1) - 15*mHu2*Sqr(g2) - 30*AbsSqr(
      MassWB)*Sqr(g2) - 10*mHu2*Sqr(gp)*Sqr(QHu) - 10*mHu2*Sqr(gp)*Sqr(Qq) - 20
      *AbsSqr(MassU)*Sqr(gp)*(Sqr(QHu) + Sqr(Qq) - Sqr(Qu)) + 10*mHu2*Sqr(gp)*
      Sqr(Qu))*(Yu*Yu.adjoint()) + 0.8*twoLoop*(MassB*Sqr(g1) - 5*(3*
      traceAdjYuTYu + traceAdjYvTYv + 3*MassWB*Sqr(g2) + 2*MassU*Sqr(gp)*(Sqr(
      QHu) + Sqr(Qq) - Sqr(Qu))) - 5*Conj(Lambdax)*TLambdax)*(Yu*(TYu).adjoint(
      )) - 0.8*twoLoop*(-(Conj(MassB)*Sqr(g1)) + 5*(3*traceconjTYuTpYu +
      traceconjTYvTpYv + Conj(TLambdax)*Lambdax + 3*Conj(MassWB)*Sqr(g2) + 2*
      Conj(MassU)*Sqr(gp)*(Sqr(QHu) + Sqr(Qq) - Sqr(Qu))))*(TYu*Yu.adjoint()) -
      0.8*twoLoop*(5*AbsSqr(Lambdax) + Sqr(g1) - 5*(-3*traceYuAdjYu -
      traceYvAdjYv + 3*Sqr(g2) + 2*Sqr(gp)*(Sqr(QHu) + Sqr(Qq) - Sqr(Qu))))*(
      TYu*(TYu).adjoint()) - 0.4*twoLoop*(5*AbsSqr(Lambdax) + Sqr(g1) - 5*(-3*
      traceYuAdjYu - traceYvAdjYv + 3*Sqr(g2) + 2*Sqr(gp)*(Sqr(QHu) + Sqr(Qq) -
      Sqr(Qu))))*(mu2*Yu*Yu.adjoint()) - 0.8*twoLoop*(5*AbsSqr(Lambdax) + Sqr(
      g1) - 5*(-3*traceYuAdjYu - traceYvAdjYv + 3*Sqr(g2) + 2*Sqr(gp)*(Sqr(QHu)
      + Sqr(Qq) - Sqr(Qu))))*(Yu*mq2*Yu.adjoint()) - 0.4*twoLoop*(5*AbsSqr(
      Lambdax) + Sqr(g1) - 5*(-3*traceYuAdjYu - traceYvAdjYv + 3*Sqr(g2) + 2*
      Sqr(gp)*(Sqr(QHu) + Sqr(Qq) - Sqr(Qu))))*(Yu*Yu.adjoint()*mu2) - 4*(mHd2
      + mHu2)*twoLoop*(Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 4*twoLoop*(Yu*
      Yd.adjoint()*TYd*(TYu).adjoint()) - 8*mHu2*twoLoop*(Yu*Yu.adjoint()*Yu*
      Yu.adjoint()) - 4*twoLoop*(Yu*Yu.adjoint()*TYu*(TYu).adjoint()) - 4*
      twoLoop*(Yu*(TYd).adjoint()*TYd*Yu.adjoint()) - 4*twoLoop*(Yu*(TYu)
      .adjoint()*TYu*Yu.adjoint()) - 4*twoLoop*(TYu*Yd.adjoint()*Yd*(TYu)
      .adjoint()) - 4*twoLoop*(TYu*Yu.adjoint()*Yu*(TYu).adjoint()) - 4*twoLoop
      *(TYu*(TYd).adjoint()*Yd*Yu.adjoint()) - 4*twoLoop*(TYu*(TYu).adjoint()*
      Yu*Yu.adjoint()) - 2*twoLoop*(mu2*Yu*Yd.adjoint()*Yd*Yu.adjoint()) - 2*
      twoLoop*(mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 4*twoLoop*(Yu*mq2*
      Yd.adjoint()*Yd*Yu.adjoint()) - 4*twoLoop*(Yu*mq2*Yu.adjoint()*Yu*
      Yu.adjoint()) - 4*twoLoop*(Yu*Yd.adjoint()*md2*Yd*Yu.adjoint()) - 4*
      twoLoop*(Yu*Yd.adjoint()*Yd*mq2*Yu.adjoint()) - 2*twoLoop*(Yu*Yd.adjoint(
      )*Yd*Yu.adjoint()*mu2) - 4*twoLoop*(Yu*Yu.adjoint()*mu2*Yu*Yu.adjoint())
      - 4*twoLoop*(Yu*Yu.adjoint()*Yu*mq2*Yu.adjoint()) - 2*twoLoop*(Yu*
      Yu.adjoint()*Yu*Yu.adjoint()*mu2) + 0.035555555555555556*twoLoop*(15*(
      -7.745966692414834*g1*(gp*Qu*(Tr2U114 + Tr2U141) + Tr31) + 15*gp*Qu*(gp*
      Qu*Tr2U144 + Tr34) + 20*Tr23*Quad(g3) + 4*Tr2U111*Sqr(g1)) + 4*AbsSqr(
      MassB)*Sqr(g1)*(321*Sqr(g1) + 80*Sqr(g3)))*UNITMATRIX(3))).real();
   const Eigen::Matrix<double,3,3> beta_mu2_2 = (0.17777777777777778*
      twoLoop*(2*Conj(MassB)*Sqr(g1)*(16*MassG*Sqr(g3) - 3*(2*MassB + MassU)*(9
      *Qd + 9*Qe - 3*QHd + 3*QHu - 9*Ql + 9*Qq - 22*Qu)*Qu*Sqr(gp)) + 8*Conj(
      MassG)*Sqr(g3)*(4*(MassB + 2*MassG)*Sqr(g1) + 15*(-2*MassG*Sqr(g3) + (2*
      MassG + MassU)*Sqr(gp)*Sqr(Qu))) + 3*Qu*Conj(MassU)*Sqr(gp)*(-2*(MassB +
      2*MassU)*(9*Qd + 9*Qe - 3*QHd + 3*QHu - 9*Ql + 9*Qq - 22*Qu)*Sqr(g1) + 5*
      Qu*(8*(MassG + 2*MassU)*Sqr(g3) + 9*MassU*Sqr(gp)*(9*Sqr(Qd) + 3*Sqr(Qe)
      + 2*Sqr(QHd) + 2*Sqr(QHu) + 6*Sqr(Ql) + 18*Sqr(Qq) + Sqr(Qs) + 11*Sqr(Qu)
      + 3*Sqr(Qv)))))*UNITMATRIX(3)).real();

   beta_mu2 = beta_mu2_1 + beta_mu2_2;


   return beta_mu2;
}

/**
 * Calculates the 3-loop beta function of mu2.
 *
 * @return 3-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_3_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return beta_mu2;
}

} // namespace flexiblesusy
