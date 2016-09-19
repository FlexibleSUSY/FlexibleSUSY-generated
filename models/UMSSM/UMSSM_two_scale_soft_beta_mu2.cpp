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

// File generated at Mon 19 Sep 2016 09:59:47

#include "UMSSM_two_scale_soft_parameters.hpp"
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
 * Calculates the one-loop beta function of mu2.
 *
 * @return one-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_one_loop(const Soft_traces& soft_traces) const
{
   const auto Qu = INPUT(Qu);
   const double Tr11 = TRACE_STRUCT.Tr11;
   const double Tr14 = TRACE_STRUCT.Tr14;


   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = (oneOver16PiSqr*(4*mHu2*(Yu*Yu.adjoint()) + 4*(TYu*(TYu)
      .adjoint()) + 2*(mu2*Yu*Yu.adjoint()) + 4*(Yu*mq2*Yu.adjoint()) + 2*(Yu*
      Yu.adjoint()*mu2) - 1.0327955589886444*g1*Tr11*UNITMATRIX(3) + 2*gp*Qu*
      Tr14*UNITMATRIX(3) - 2.1333333333333333*AbsSqr(MassB)*Sqr(g1)*UNITMATRIX(
      3) - 10.666666666666666*AbsSqr(MassG)*Sqr(g3)*UNITMATRIX(3) - 8*AbsSqr(
      MassU)*Sqr(gp)*Sqr(Qu)*UNITMATRIX(3))).real();


   return beta_mu2;
}

/**
 * Calculates the two-loop beta function of mu2.
 *
 * @return two-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_two_loop(const Soft_traces& soft_traces) const
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
      )) - 0.008888888888888889*twoLoop*(1350*traceconjTYuTpYu + 450*
      traceconjTYvTpYv + 450*Conj(TLambdax)*Lambdax - 90*Conj(MassB)*Sqr(g1) +
      1350*Conj(MassWB)*Sqr(g2) + 900*Conj(MassU)*Sqr(gp)*Sqr(QHu) + 900*Conj(
      MassU)*Sqr(gp)*Sqr(Qq) - 900*Conj(MassU)*Sqr(gp)*Sqr(Qu))*(TYu*Yu.adjoint
      ()) - 0.008888888888888889*twoLoop*(1350*traceYuAdjYu + 450*traceYvAdjYv
      + 450*AbsSqr(Lambdax) + 90*Sqr(g1) - 1350*Sqr(g2) - 900*Sqr(gp)*Sqr(QHu)
      - 900*Sqr(gp)*Sqr(Qq) + 900*Sqr(gp)*Sqr(Qu))*(TYu*(TYu).adjoint()) -
      0.008888888888888889*twoLoop*(675*traceYuAdjYu + 225*traceYvAdjYv + 225*
      AbsSqr(Lambdax) + 45*Sqr(g1) - 675*Sqr(g2) - 450*Sqr(gp)*Sqr(QHu) - 450*
      Sqr(gp)*Sqr(Qq) + 450*Sqr(gp)*Sqr(Qu))*(mu2*Yu*Yu.adjoint()) -
      0.008888888888888889*twoLoop*(1350*traceYuAdjYu + 450*traceYvAdjYv + 450*
      AbsSqr(Lambdax) + 90*Sqr(g1) - 1350*Sqr(g2) - 900*Sqr(gp)*Sqr(QHu) - 900*
      Sqr(gp)*Sqr(Qq) + 900*Sqr(gp)*Sqr(Qu))*(Yu*mq2*Yu.adjoint()) -
      0.008888888888888889*twoLoop*(675*traceYuAdjYu + 225*traceYvAdjYv + 225*
      AbsSqr(Lambdax) + 45*Sqr(g1) - 675*Sqr(g2) - 450*Sqr(gp)*Sqr(QHu) - 450*
      Sqr(gp)*Sqr(Qq) + 450*Sqr(gp)*Sqr(Qu))*(Yu*Yu.adjoint()*mu2) -
      0.008888888888888889*(450*mHd2 + 450*mHu2)*twoLoop*(Yu*Yd.adjoint()*Yd*
      Yu.adjoint()) - 4*twoLoop*(Yu*Yd.adjoint()*TYd*(TYu).adjoint()) - 8*mHu2*
      twoLoop*(Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 4*twoLoop*(Yu*Yu.adjoint()*
      TYu*(TYu).adjoint()) - 4*twoLoop*(Yu*(TYd).adjoint()*TYd*Yu.adjoint()) -
      4*twoLoop*(Yu*(TYu).adjoint()*TYu*Yu.adjoint()) - 4*twoLoop*(TYu*
      Yd.adjoint()*Yd*(TYu).adjoint()) - 4*twoLoop*(TYu*Yu.adjoint()*Yu*(TYu)
      .adjoint()) - 4*twoLoop*(TYu*(TYd).adjoint()*Yd*Yu.adjoint()) - 4*twoLoop
      *(TYu*(TYu).adjoint()*Yu*Yu.adjoint()) - 2*twoLoop*(mu2*Yu*Yd.adjoint()*
      Yd*Yu.adjoint()) - 2*twoLoop*(mu2*Yu*Yu.adjoint()*Yu*Yu.adjoint()) - 4*
      twoLoop*(Yu*mq2*Yd.adjoint()*Yd*Yu.adjoint()) - 4*twoLoop*(Yu*mq2*
      Yu.adjoint()*Yu*Yu.adjoint()) - 4*twoLoop*(Yu*Yd.adjoint()*md2*Yd*
      Yu.adjoint()) - 4*twoLoop*(Yu*Yd.adjoint()*Yd*mq2*Yu.adjoint()) - 2*
      twoLoop*(Yu*Yd.adjoint()*Yd*Yu.adjoint()*mu2) - 4*twoLoop*(Yu*Yu.adjoint(
      )*mu2*Yu*Yu.adjoint()) - 4*twoLoop*(Yu*Yu.adjoint()*Yu*mq2*Yu.adjoint())
      - 2*twoLoop*(Yu*Yu.adjoint()*Yu*Yu.adjoint()*mu2) - 0.008888888888888889*
      twoLoop*(-1200*Power(g3,4)*Tr23*UNITMATRIX(3) + 464.75800154489*g1*gp*Qu*
      Tr2U114*UNITMATRIX(3) + 464.75800154489*g1*gp*Qu*Tr2U141*UNITMATRIX(3) +
      464.75800154489*g1*Tr31*UNITMATRIX(3) - 900*gp*Qu*Tr34*UNITMATRIX(3) -
      5136*Power(g1,4)*AbsSqr(MassB)*UNITMATRIX(3) - 240*Tr2U111*Sqr(g1)*
      UNITMATRIX(3) - 1280*AbsSqr(MassB)*Sqr(g1)*Sqr(g3)*UNITMATRIX(3) - 900*
      Tr2U144*Sqr(gp)*Sqr(Qu)*UNITMATRIX(3)))).real();
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
 * Calculates the three-loop beta function of mu2.
 *
 * @return three-loop beta function
 */
Eigen::Matrix<double,3,3> UMSSM_soft_parameters::calc_beta_mu2_three_loop(const Soft_traces& soft_traces) const
{
   DEFINE_PROJECTOR(3,3,3,3)



   Eigen::Matrix<double,3,3> beta_mu2;

   beta_mu2 = ZEROMATRIX(3,3);


   return beta_mu2;
}

} // namespace flexiblesusy
